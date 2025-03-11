import numpy as np

def LDPC_decoder(c_llr, H, iter_num):
  """
    Implements the classical sum-product algorithm (SPA) for LDPC codes.

    Parameters:
    - c_llr (numpy.ndarray): Channel LLR (Log-Likelihood Ratio) values.
    - H (numpy.ndarray): Parity check matrix.
    - iter_num (int): Total number of iterations for decoding.

    Returns:
    - m_est (numpy.ndarray): Estimated message bits after decoding.

    This function performs LDPC decoding using the Sum-Product Algorithm (SPA). It takes the channel LLR values,
    parity check matrix, and the number of iterations as inputs, and returns the estimated message bits after decoding.

    Example:
    ```
    # Define inputs
    c_llr = np.array([0.1, -0.3, 0.5, -0.2])
    H = np.array([[1, 0, 1, 0],
                  [0, 1, 1, 1]])
    iter_num = 10

    # Perform LDPC decoding
    m_est = LDPC_decoder(c_llr, H, iter_num)
    ```
    """

  size_H = H.shape
  m,n = size_H
  
  # Initiliaze the messages - 
  # For the first initial message passing,
  # realize with the channel values.
  VN_to_CN = H*c_llr 

  for iter in range(1,iter_num+1):
      CN_to_VN = np.zeros(size_H)
    
      # Calculate CN to VN messages
      for i in range(m):

          # Get current row vector of H
          curr_row = H[i,:]
          # Find where this equals to non-zero
          ind = np.squeeze(np.argwhere(curr_row))
          # get messages coming to CN
          mes = VN_to_CN[i,ind]

          for j in range(len(ind)):

              # Delete the corresponding CN_to_VN meesage where we send info
              mes_upd = np.copy(mes)
              mes_upd = np.delete(mes_upd, j)

              # Find current CN to VN message
              CN_to_VN[i, ind[j]] = 2*np.arctanh(np.clip(np.prod(np.tanh(mes_upd/2)), -0.9999999999, 0.9999999999))

      # Estimate the messages
      llr = np.sum(CN_to_VN, 0) + c_llr
      m_est = llr < 0
      if np.sum(np.mod(m_est @ H.T, 2)) == 0:
          #print('Process is stopped at iteration: ', iter)
          break;

      VN_to_CN = np.zeros(size_H)
      # Calculate VN to CN messages
      for j in range(n):

          # Get current column vector of H
          curr_col = H[:,j]
          # Find where this equals to non-zero
          ind = np.squeeze(np.argwhere(curr_col))
          # get messages coming to CN
          mes = CN_to_VN[ind, j]

          for i in range(len(ind)):

              # Delete the corresponding CN_to_VN meesage where we send info
              mes_upd = np.copy(mes)
              mes_upd = np.delete(mes_upd, i)

              # Find current VN to CN message
              VN_to_CN[ind[i],j] = np.sum(mes_upd)+c_llr[0,j]

  # Estimate the messages
  llr = np.sum(CN_to_VN, 0) + c_llr
  m_est = llr < 0

  return m_est, llr