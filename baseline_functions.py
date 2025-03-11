import numpy as np
from insertion_deletion import ins_del_channel, insert_regular_markers
import itertools
from reciever import Receiver_Regular
import sys

def baseline_simulation_LDPC(H, iter_num, marker_seq, Nc, Pd, Pi, Ps, 
                            print_every = 10, max_fer = 100, max_codes = 10000,
                            CSI='known', P_est=None, filepath = None):

    """
    Perform baseline decoder testing for LDPC decoding using the BCJR algorithm 
    with insertion/deletion/substitution channels.

    Args:
        H (numpy.ndarray): Parity check matrix of the LDPC code.
        iter_num (int): Number of iterations for LDPC decoding (sum-product algorithm).
        marker_seq (array_like): Marker sequence.
        Nc (int): Total codeword block size.
        Pi (list): Insertion test probability of the channel(s).
        Ps (list): Substitution test probability of the channel(s).
        Pd (list): Deletion test probability of the channel(s).
        print_every (int): Frequency of printing progress during testing.
        CSI(str, choices = {'known', 'unk}): Channel state information (channel probabilities) is known or not
        P_est(list): Channel probability estimates if CSI is set to 'unk'. [Pi_est, Ps_est]

    Returns:
        res

    Note:
    This function performs baseline testing for LDPC decoding by simulating transmission over
    insertion/deletion/substitution channels. It generates random messages, encodes them using
    regular markers, applies the specified channel model (with insertion, deletion and substition probs.), 
    performs LDPC decoding after estimating LLRs with the BCJR algorithm, 
    and calculates the error rates (BER, FER).
    """

    # Check if correct inputs are taken !
    if  CSI != 'unk' and P_est == None:
        assert CSI != 'unk' and P_est == None, 'If CSI is unknown (unk), then P_est should not be empty!'
    elif CSI != 'unk' and CSI != 'known':
        assert CSI != 'known' and CSI != 'unk', "CSI should be set to 'known' or 'unk'!"

    # Codeword Lenght
    n = H.shape[-1]

    # Number of padded bits
    num_padded_bits = Nc - (n%Nc)

    # Create interleaver (permutation)
    if n % Nc != 0:
        perm = np.random.permutation(n + num_padded_bits)
    else:
        perm = np.random.permutation(n)

    # Total codeword block
    Nr = marker_seq.shape[-1]
    N = Nc + Nr
    
    # Create testing points
    test_points = list(itertools.product(Pd,Pi,Ps))

    print(f"------------SIMULATION STARTS---------------")
    
    # Hold results array.
    res = []
    for i, test_point in enumerate(test_points):

        # Init BER/FER error metrics
        ber_total = 0
        fer_total = 0

        # Get current Pd, Pi, Ps
        Pd_point, Pi_point, Ps_point = test_point
        code_simulated = 1
        print(f"------------Test point (Pd, Pi, Ps) = {test_point}---------------")
        while fer_total <= max_fer and code_simulated <= max_codes:

            # Create a sequence and add it to all zero codeword.
            # This step is done in order to bypass the creation of the codeword via
            # generator matrix of the LDPC code which is harder to obtain for some
            # cases.
            rand_seq = np.random.randint(low = 0, high=2, size=(1,n), dtype=int)
            m = np.zeros((1,n)) + rand_seq

            # Insert markers
            c,mask = insert_regular_markers(m, Nc, marker_seq, perm)
            T = c.shape[-1]

            # Send codeword through the general insertion/deletion/substition channel.
            y, _ = ins_del_channel(c, Pd_point, Pi_point, Ps_point)
            R = len(y)

            if CSI == 'known': # If channel probabilities are known
                m_est = Receiver_Regular(y, T, Pd_point, Pi_point, Ps_point, mask, marker_seq, 
                                                rand_seq, iter_num, n, perm, N, Nc, H)
            elif CSI == 'unk': # If channel probabilities are not known
                Pd_est = (T-R)/T

                if Pd_est == 0: 
                    Pd_est = sys.float_info.min

                m_est = Receiver_Regular(y, T, Pd_est, P_est[0],P_est[1], mask, marker_seq,  
                                                rand_seq, iter_num, n, perm, N, Nc, H)


            # Add errors
            ber = np.sum(m_est != np.zeros((1,n)))
            ber_total += ber

            if ber > 0:
                fer_total += 1

            if code_simulated % print_every == 0:
                ber = ber_total/(n*code_simulated)
                fer = fer_total/code_simulated
                print(f"{code_simulated}) (Pd, Pi, Ps) = {test_point}, BER: {ber: .6f}, FER: {fer: .6f}")

            code_simulated += 1

        print(f"{i+1}) (Pd, Pi, Ps) = {test_point} testing finished, BER: {ber_total/(n*code_simulated) : .6f}, FER: {fer_total/code_simulated : .6f}")
        
        # Append results.
        res.append((ber_total/(code_simulated*n),fer_total/(code_simulated)))
    
    # Save the results.
    return res

