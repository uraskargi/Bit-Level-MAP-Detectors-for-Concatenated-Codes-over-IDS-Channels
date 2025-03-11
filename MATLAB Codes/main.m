% Clear all outputs & start
clc; clear; close all; format long;

% Load Parity Check Matrix -
% Obtained from David Mackay's page
% (https://www.inference.org.uk/mackay/CodesFiles.html)
% -------------------------------------------------------------
% We are using the spesific regular LDPC code used in Prof. Duman's paper
% link here for that code ----> https://www.inference.org.uk/mackay/codes/data.html#l141
% n = 16383, k = 14252, m = 2131, column weight(t) = 4, row weight = 31.
% If you want to use this matrix, please call
% load('H_matrix_large.mat')
% -------------------------------------------------------------
% If you want to another (Regular LDPC code, n = 804, k = 402)
% call load('H_matrix_medium.mat'). 
% With this smaller code, we can run the code
% much faster than the first one. 
% -------------------------------------------------------------
% If you want to another (Regular LDPC code, n = 204, k = 102)
% call load('H_matrix_small.mat'). 
% With this smaller code, we can run the code
% much faster than the first one. 
% -------------------------------------------------------------
% If you want to another (n = 32000, m = 2240)
% call load('H_matrix_xlarge.mat'). 
% With this smaller code, we can run the code
% much faster than the first one. 
% -------------------------------------------------------------

% ALL PARAMETERS OF THE CODE
load('Matrices/H_new.mat')    % Load the parity-check matrix for LDPC code
iter_LDPC = 100;                         % Total iteration number for the LDPC 
                                            % SPA (sum-product algoritm) based decoder
early_finish = 1;                       % Early finish for LDPC decoder. It means that
                                            % when llr @ VN produces a
                                            % valid codeword, then decoder
                                            % is finished.              
iter_block = 1;                         % Total iteration number for the block
seed = 100;                             % Seed for reproducibility
Nc = 12;                                % Specify Nc
marker_code = [0,1,0];                    % Specify the marker sequence in this array. 
                                            % The same marker consisting of 𝑁m 
                                            % consecutive bits is inserted 
                                            % every 𝑁c bits at the output of t
                                            % he outer encoder.
Pi_test = [0.005];      % Define testing insertion probs - in an array
Pd_test = [0.01,0.02,0.03,0.04,0.05];           % Define testing deletion probs -  in an array 
Ps_test = [0.005];               % Define testing subs probs - in an array
print_every = 10;                      % Show results for every print_every codewords

max_fer = 200;                         % The frame errors total when simukation ends!
                                            % Specify so that simulation
                                            % ends when Frame errors reach
                                            % max_fer
max_codes = 50000;                     % The frame errors total when simukation ends!
                                            % Specify so that simulation
                                            % ends when total number of
                                            % codewords simulated reach
                                            % max_codes
start_with_marker = 0;

estimation = 0;                        % Esimation = 0 means that channel probabilities
                                            % are perfectly known to the
                                            % receiver. Estimation = 1
                                            % means that channel
                                            % probabilities are not known
                                            % at the recevier
Ps_est = [];
Pi_est = [];
Pd_est = [];
% Set the seed for reproducibility.
rng(seed)

% These lines below create all testing points for the simulation.
[d, i, s] = ndgrid(Pd_test, Pi_test, Ps_test);
test_points = [d(:) i(:) s(:)];
num_test_points = size(test_points,1);


% Useful variables !
size_H = size(H);                       % Size of H
n = size_H(2);                          % n: Codeword length
m = size_H(1);                          % m: Number of check nodes,i.e, m=n-k
k = n-m;                                % k: Message length
R = k/n;                                % rate(of outer code)= k/n
Nm = length(marker_code);               % Nm = Total number of marker bits
                                        %       in every Nc bits.
N = Nc + Nm;                            % Total marker and meesage bit block
rm = Nc/(Nc+Nm);                        % Rate of marker code (inner code)
mtotal = k;                             % Total message bits for each block
num_padded_bits = Nc - rem(n, Nc);      % Number of padded bits

% Create random permuation for interleaver
% if total codeword lenght does not divide Nc, we add some padding 0's
if  rem(n, Nc) ~= 0 
    %perm = randperm(n + num_padded_bits); 
    perm = 1:n + num_padded_bits;
else
    %perm = randperm(n); 
    perm = 1:n;
end

res = [];   % This array holds the BER & FER results 

disp('---------Simulation Starts---------')
for z = 1:num_test_points

    % Init. bit error metric
    bit_errors = 0;
    frame_errors = 0;
    num_code = 0;
       
    % Get current Simulation parameters
    Pd = test_points(z,1);
    Pi = test_points(z,2);
    Ps = test_points(z,3);
    
    % Simulation starts for the test point
    info = [
            sprintf('Simulation starts for the test point (Pd=%.3f,', Pd)...
            sprintf('Pi=%.3f,', Pi)...
            sprintf('Ps=%.3f)', Ps)];
    disp(info)

    while (frame_errors < max_fer) && (num_code < max_codes)
        
        % A-priori LLR !!!
        ext_llr = zeros(length(perm), 1);
         
        % I wrote this line to correct someeeeeeeettttthhhiiiiings :=)
        % Create codeword
        rand_seq = randi([0,1],1,n);
        m = rand_seq;
        
        % Create marker coded bit and recieve codeword and mask vector
        % mask vector --> 1 equals when marker bit indices are observed
        [c, mask]  = create_ldpc_marker_code(m, Nc, marker_code, perm, start_with_marker);
        % Send it through the insertion deletion channel
        [y, ~] = ins_del_channel(c, Pd, Pi, Ps);
        
        % lengths of TRANSMITTED vecors etc
        T = length(c);
        R = length(y);
        
        % This part is included for CSI
        if estimation == 1
            Pd_alg = Pd_est;
            Ps_alg = Ps_est;
            Pi_alg = Pi_est;
        else 
            Pd_alg = Pd;
            Ps_alg = Ps;
            Pi_alg = Pi;
        end
        
        for f = 1:iter_block

            % Forward Recursion
            alpha = forward_recursion_log_domain(y,ext_llr,T,Pd_alg, Pi_alg,...
                Ps_alg,mask,marker_code,N,Nc);
            % Backward Recursion
            beta = backward_recursion_log_domain(y,ext_llr,T,Pd_alg, Pi_alg,...
                Ps_alg,mask,marker_code,N,Nc);
            % Message from forward & backward Equations
            [c_est, llr] = estimate_message(y, alpha, beta,Pd_alg,Pi_alg, ...
                Ps_alg, T);

            llr = llr(mask == 0);              % Get rid of the llrs of the non-marker bits                             
            llr2 = zeros(size(llr));
            llr2(perm) = llr;                  % Apply de-interleaver 
            llr2 = llr2(1:n);                  % Get rid of the padded bits
            llr2(rand_seq == 1) = -llr2(rand_seq == 1);  % Change the effect of the random sequence
            
            
            % Estimate the message bits (or get LLRs) from lDPC decoder                       
            [m_est, ext_llr] = LDPC_decoder(llr2, H, iter_LDPC, 0);
    
            % Pad zero LLRs for randomly added bits to the orginal message
            ext_llr = [ext_llr zeros(1,num_padded_bits)];     % Add 0 llr's for padded bits ! 
            ext_llr(rand_seq == 1) = -ext_llr(rand_seq == 1); % Change the sign for 1's
            ext_llr = ext_llr(perm) ;                         % Permute (Interleve)
                                
        end
        
        % Update BER and FER metrics
        ber = sum(m_est ~= zeros(1,n));
        bit_errors = bit_errors + ber;
        frame_errors = frame_errors + (1-isequal(zeros(1,n),m_est));
        num_code = num_code + 1;   % Increae Codeword Number
        
        % Print the results once every
        if rem(num_code,print_every) == 0
            results = [
            sprintf('N: %d - ' , num_code) ...
            sprintf('(Pd=%.3f,', Pd)...
            sprintf('Pi=%.3f,', Pi)...
            sprintf('Ps=%.3f) ', Ps)...
            sprintf('BER=%.6f ', bit_errors/(n*num_code)) ...
            sprintf('FER=%.6f ', frame_errors/num_code)];
            disp(results)
        end
    end

% Calculate bit error rate
ber = bit_errors/(num_code*n);
fer = frame_errors/num_code;
res = [res [ber;fer]];
results = [
            sprintf('Simulation Ends at the iteration: %d - ' , num_code) ...
            sprintf('(Pd=%.3f,', Pd)...
            sprintf('Pi=%.3f,', Pi)...
            sprintf('Ps=%.3f) ', Ps)...
            sprintf('BER=%.6f ', bit_errors/(n*num_code)) ...
            sprintf('FER=%.6f ', frame_errors/num_code)];
 disp(results)
end





