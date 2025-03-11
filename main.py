import os
import numpy as np
import argparse
import scipy
import sys
from baseline_functions import baseline_simulation_LDPC
from utils import save_args_to_file

parser = argparse.ArgumentParser(description='Testing for the baseline decoders!')

parser.add_argument('--main-dir', type=str, default= "./Results")
parser.add_argument('--test_points_Pd', default=[0.04],
                    type=list, help = 'Deletion probability test points (type: list)')
parser.add_argument('--test_points_Ps', default=[sys.float_info.min], type=list, 
                    help = 'Substitution probability test points (type: list)')
parser.add_argument('--test_points_Pi', default=[sys.float_info.min], type=list, 
                    help = 'Insertion probability test points (type: list)')
parser.add_argument('--fer_errors', default=100, type=int, 
                    help = 'Maximum number of FER errors for each test point.')       
parser.add_argument('--marker_sequence', default=np.array([0,1,0]).reshape(1,-1), type = np.array,
                      help = 'Specify the marker sequence (with shape format [1,-1]).')   
parser.add_argument('--Nc', default=10, type = int, help = 'After how many bits, markers are inserted')
parser.add_argument('--iter_num', default=100, type = int, help='Iteration number of the LDPC sum product Decoder')
parser.add_argument('--matrix', choices = ['small', 'medium', 'large', 'xlarge'], type = str,
                    help = 'The size of the parity check matrix (default small). Experiments show results only for small', 
                    default = 'small')
parser.add_argument('--print_every', default=10, type = int)
parser.add_argument('--CSI', choices=['unk', 'known'], help = 'Channel State Information', default = 'unk')
parser.add_argument('--seed', default=1000, type = int, help = 'Seed for reproducibility')

def main():

    # Fetch the arguments!
    args = parser.parse_args()
    # Fix the seed for reproducibility!
    np.random.seed(args.seed)

    # Fix the seed for reproducibility!
    matrix_path = 'H_matrix_' + args.matrix + '.mat'
    filepath = os.path.join('Matrices', matrix_path)
    H = scipy.io.loadmat(filepath)['H'] 
    main_dir = args.main_dir
    k = H.shape[0]
    n = H.shape[0]

    # Get args
    iter_num = args.iter_num
    marker_seq = args.marker_sequence
    print_every = args.print_every
    Nc = args.Nc
    Nm = marker_seq.shape[-1]
    fer_errors = args.fer_errors
    Pi = args.test_points_Pi
    Ps = args.test_points_Ps
    Pd = args.test_points_Pd
    #path = args.path
    csi = args.CSI
    #csi = 'unk'

    # Simulation logs
    print('Simulation starts for the baseline decoder with BCJR based estimator and the LDPC decoder')
    print('The Partiy-Check matrix size is ', H.shape)
    print('Nc (In how many codeword bits, marker sequence is added) =', Nc)
    print('Iteration Number of the LDPC decoder = ', iter_num)
    print('Marker sequence is ', marker_seq)

    save_dir = "Nc_%d_Nm_%d_k_%d_n_%d" % (Nc, Nm,k,n)
    # Directories
    main_dir2 = os.path.join(main_dir, args.matrix)
    directory = os.path.join(main_dir2, save_dir)

    if not os.path.exists(main_dir):
        os.mkdir(main_dir)
    if not os.path.exists(main_dir2):
        os.mkdir(main_dir2)
    if not os.path.exists(directory):
        os.mkdir(directory)

    # Save parameters
    save_args_to_file(args, directory)

    # Run the simulation
    res = baseline_simulation_LDPC(H,iter_num,marker_seq,Nc,Pd,Pi,Ps,
                                   print_every = print_every,
                                    CSI=csi, P_est=[sys.float_info.min, sys.float_info.min], 
                                    max_fer = fer_errors)
    
    test_directory = os.path.join(directory, 'test_results')
    stri = "test_results"
    if not os.path.exists(test_directory):
            os.mkdir(test_directory )
    np.save(os.path.join(test_directory, stri), res)
    print('Simulation is finished')


if __name__ == '__main__':
  main()

