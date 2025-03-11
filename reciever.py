import numpy as np
from BCJR import BCJR_Regular
from decoders import LDPC_decoder
 
def Receiver_Regular(y, T, Pd, Pi, Ps, mask, marker_seq, rand_seq, iter_num, n, perm, N, Nc, H):
    
    # BCJR
    llr = BCJR_Regular(y,T,Pd,Pi,Ps,mask,marker_seq,N,Nc)

    # Get rid of the llrs corresponding to marker bits
    llr_new = llr[np.squeeze(mask) == 0, :].T

    # Deinterleve
    llr = np.zeros(llr_new.shape)
    llr[:,perm] = llr_new

    # Get rid of the extra bits padded at the end of seqeunce while creating codewords
    llr = llr[:,0:n]
    
    # This step is done in order to reverse the effect of the random sequence added 
    # at the beginning.
    llr[:,np.squeeze(rand_seq) == 1] = -llr[:,np.squeeze(rand_seq) == 1]

    # LDPC Decoder (clip values in order to avoid numerical errors) 
    llr = np.clip(llr, -30, 30)
    m_est, _ = LDPC_decoder(llr,H,iter_num)

    return m_est