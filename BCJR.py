import numpy as np
from utils import max_star

def BCJR_Regular(y,T,Pd,Pi,Ps,mask,marker_seq,N,Nc): 

    # Calculate forward probabilities.
    alpha = forward_recursion_log_domain(np.squeeze(y), T, Pd, Pi, Ps, np.squeeze(mask), 
                                        np.squeeze(marker_seq), N, Nc)
    # Calculate backward probabilities.                                 
    beta = backward_recursion_log_domain(np.squeeze(y), T, Pd, Pi, Ps, np.squeeze(mask), 
                                        np.squeeze(marker_seq), N, Nc)
    # Estimate LLRs from forward and backward probabilities.
    _, llr = estimate_message(np.squeeze(y), alpha, beta, Pd, Pi, Ps, T)

    return llr

def forward_recursion_log_domain(y, T, Pd, Pi, Ps, mask, marker_seq, N, Nc):
    """
    Calculates forward probabilities in BCJR algorithm deployed for 
    insertion/deletion/substitution channels.
    
    Args:
        y (list): Channel realization (Output of the channel).
        T (int): Total transmitted bits.
        Pd (float): Deletion probability of the Channel.
        Pi (float): Insertion probability of the Channel.
        Ps (float): Substitution probability of the Channel.
        mask (array_like): Indicates whether the current bit is a marker bit (1 if the current bit is a marker bit.
        marker_seq (array_like): Code of the marker.
        N (int): Total block length (N = Nc + Nm).
        Nc (int): Total codeword block.

    Returns:
        alpha(np.darray): Forward probabilities.
    """

    R = len(y)                              # Received bit numbers 
    alpha = -np.inf * np.ones((R+1, T+1))   # Forward recursion inits
    alpha[0, 0] = 0                         # Init. the starting point
    Pt = 1 - Pd - Pi                        # Transmission Prob.

    for k in range(1, T+1): # From stage 1 to T
        
        for n in range(max(0, R-2*(T-k)), min(2*k+1, R+1)):
            
            if n == 0: # For state n = 0, only deletion event occurs
                alpha[n, k] = np.log(Pd) + alpha[n, k-1]
                
            elif n == 1: # For state n = 1, deletion or transmission events can occur
                a = np.log(Pd) + alpha[n, k-1]
                # If this is a marker bit position, we need the marker bit @ that position
                if mask[k-1] == 1:
                    curr_marker_bit = marker_seq[(k-1)%N - Nc]
                    b = np.log(Pt) + np.log((1-2*Ps)*(curr_marker_bit == y[n-1]) + Ps) + alpha[n-1, k-1]
                else:
                    b = np.log(Pt/2) + alpha[n-1, k-1]
                alpha[n, k] = max_star(a, b)

            elif n > 1: # For states bigger than 1, deletion, transmission, or insertion can occur
                a = np.log(Pd) + alpha[n, k-1]
                b = np.log(Pi/4) + alpha[n-2, k-1]
                if mask[k-1] == 1:
                    curr_marker_bit = marker_seq[(k-1)%N - Nc]
                    c = np.log(Pt) + np.log((1-2*Ps)*(curr_marker_bit == y[n-1]) + Ps) + alpha[n-1, k-1]
                else:
                    c = np.log(Pt/2) + alpha[n-1, k-1]
                alpha[n, k] = max_star(a, b, c)
    return alpha


def backward_recursion_log_domain(y, T, Pd, Pi, Ps, mask, marker_seq, N, Nc):
    """
    Calculates backward probabilities in BCJR algorithm deployed for 
    insertion/deletion/substitution channels.
    
    Args:
        y (list): Channel realization (Output of the channel).
        T (int): Total transmitted bits.
        Pd (float): Deletion probability of the Channel.
        Pi (float): Insertion probability of the Channel.
        Ps (float): Substitution probability of the Channel.
        mask (array_like): Indicates whether the current bit is a marker bit (1 if the current bit is a marker bit.
        marker_seq (array_like): Code of the marker.
        N (int): Total block length (N = Nc + Nm).
        Nc (int): Total codeword block.

    Returns:
       beta(): Backward probabilities.
    """
    
    R = len(y)                            # Received bit number
    beta = -np.inf * np.ones((R+1, T+1))  # backward recursion inits
    beta[R, T] = 0                        # Init. the starting point
    Pt = 1 - Pd - Pi

    for k in range(T, 0, -1):
        k_star = T + 1 - k
        
        # For states between R-2 to 1
        for n in range(R, max(R-1-2*k_star, -1), -1):
            
            # for states less than R - 1
            if n < R - 1:
                a = np.log(Pd)   + beta[n, k]
                b = np.log(Pi/4) + beta[n+2, k]

                if mask[k-1] == 1:
                    curr_marker_bit = marker_seq[(k-1) % N - Nc]
                    c = np.log(Pt) + np.log((1-2*Ps)*(curr_marker_bit == y[n]) + Ps) + beta[n+1, k]
                else:
                    c = np.log(Pt/2) + beta[n+1, k]

                beta[n, k-1] = max_star(a, b, c)
            
            # For state n = R - 1
            elif n == R-1:
                a = np.log(Pd) + beta[R - 1, k]
                if mask[k-1] == 1:
                    curr_marker_bit = marker_seq[(k-1) % N - Nc]
                    b = np.log(Pt) + np.log((1-2*Ps)*(curr_marker_bit == y[n]) + Ps) + beta[R, k]
                else:
                    b = np.log(Pt/2) + beta[R, k]
                beta[R-1, k-1] = max_star(a, b)

            # For state n = R
            elif n == R:
                beta[n, k-1] = np.log(Pd) + beta[n, k]

    return beta


def estimate_message(y, alpha, beta, Pd, Pi, Ps, T):
    """
    Estimate the transmitted message based on the BCJR algorithm by using
    forward and backward probabilities calculated based on the recevied sequence
    and the marker bit locations.

    Args:
        y (list): Channel realization (output of the channel).
        alpha (numpy.ndarray): Forward probabilities.
        beta (numpy.ndarray): Backward probabilities.
        Pd (float): Deletion probability of the channel.
        Pi (float): Insertion probability of the channel.
        Ps (float): Substitution probability of the channel.
        T (int): Total transmitted bits.

    Returns:
        tuple: A tuple containing:
            - m_est (numpy.ndarray): Estimated transmitted message.
            - llr (numpy.ndarray): Log-likelihood ratio of the estimated message.

    Note:
        This function estimates the transmitted message based on the forward and backward probabilities
        calculated using the BCJR algorithm. It computes the log-likelihood ratio (LLR) for each bit and
        returns the estimated message and LLR.

    """

    # Total number of received bits
    R = len(y) 

    # Transmission prob.
    Pt = 1-Pd-Pi
    
    # Likelihood arrays for individual bits
    l0 = np.zeros((T,1))
    l1 = np.zeros((T,1))

    for k in range(T):
        # Get alpha vector corresponding to alpha vector @ stage k
        alpha_k = alpha[:, k]
        # Get beta vector corresponding to beta vector @ stage (k+1)
        beta_k = beta[:, k+1]
        # Init the ins, del, trans0 and trans1 vectors to minus inf
        # Instead of minus inf, start it with a very small negative integer
        ins =     -1e100
        del_val = -1e100
        trans0 =  -1e100 
        trans1 =  -1e100 
        
        for n in range(1, min(2*k+2, R+1)):
            
            if n == R + 1:
                del_val = max_star(del_val, np.log(Pd) + beta_k[n-1] + alpha_k[n-1])
                
            elif n == R:
                del_val = max_star(del_val, np.log(Pd) + beta_k[n-1] + alpha_k[n-1])
                trans0 = max_star(trans0, np.log(Pt) + beta_k[n] + alpha_k[n-1] +
                             np.log((1-Ps)*(y[n-1] == 0) + Ps))
                trans1 = max_star(trans1, np.log(Pt) + beta_k[n] + alpha_k[n-1] +
                             np.log((1-Ps)*(y[n-1] == 1) + Ps))
            else:
                ins = max_star(ins, np.log(Pi/4) + beta_k[n+1] + alpha_k[n-1])
                del_val = max_star(del_val, np.log(Pd) + beta_k[n-1] + alpha_k[n-1])
                trans0 = max_star(trans0, np.log(Pt) + beta_k[n] + alpha_k[n-1] +
                             np.log((1-Ps)*(y[n-1] == 0) + Ps))
                trans1 = max_star(trans1, np.log(Pt) + beta_k[n] + alpha_k[n-1] +
                             np.log((1-Ps)*(y[n-1] == 1) + Ps))
                
        l0[k,0] = max_star(trans0, del_val, ins)
        l1[k,0] = max_star(trans1, del_val, ins)

    # Calculate LLR for 1,2,...,T
    llr = l0 - l1
    
    # Message estimates (It is provided in addition to LLR.)
    m_est = llr < 0

    return m_est, llr