
function beta = backward_recursion_log_domain(y, llr_apriori, T, Pd, Pi, Ps, mask, marker_code, N, Nc)
%BACKWARD_RECURSION_LOG_DOMAIN Calculates backward probabilities for
%   bit-level marker code BCJR coding in the log domain.
%
%   Inputs:
%       y           - Channel realization
%       llr_apriori - A priori log-likelihood ratios (LLRs) about information bits x_i, i = 1, 2, ..., n
%       T           - Total number of transmitted bits
%       Pd          - Deletion probability
%       Pi          - Insertion probability
%       Ps          - Substitution probability
%       mask        - Indicator vector where 1 indicates the current bit is a marker bit
%       marker_code - Code sequence of the marker
%       N           - Total block length (N = Nc + Nm)
%       Nc          - Total codeword block length
%
%   Outputs:
%       beta        - Backward probabilities in the log domain
%
%   This function implements the backward recursion in the log domain for
%   bit-level marker code BCJR coding. It calculates the backward
%   probabilities (beta) for each received bit position and transmitted
%   bit position, taking into account the a priori information, deletion,
%   insertion, and substitution events.

% Get a-priori probs from LLR information
Px_1 = flip(1./(1+exp(llr_apriori)));       % P(x_i = 1) for i = 1,2,...,n


R = length(y);                              % Recieved bit number
beta = -realmax*ones(R+1,T+1);              % Backward recursion inits
beta(R+1,T+1) = 0;                          % Init. the starting point
ind = 1;                                    % Denotes the index of the message bit
Pt = 1-Pd-Pi;                               % Transmission probability
ind_marker =  1;
flip_marker_code = flip(marker_code);
Nr = N-Nc;

for k = T:-1:1 
    k_star = T+1 - k;
     % THIS PART IS NEW - ADDED SINCE WE NEED TO TRACK THE INDEX OF THE 
    % LLR OF A-PRIORI INFORMATION. WE NEED SINCE LLR ARRAY CONTAINS INFO
    % ABOUT ONLY THE TRANSMITTED MESSAGE BITS
    if mask(k) == 0
        curr_Px_1 = Px_1(ind);
        ind = ind + 1;
    else
        if ind_marker == Nr + 1
            ind_marker = 1;
        end
        curr_marker_bit = flip_marker_code(ind_marker);  
        ind_marker = ind_marker + 1;
    end

    
    % For states betweeen R-1 to 1
    for n = min(R-1,2*T):-1:max(R+1-2*k_star,1)
        
        if n < R
            a = log(Pd)   + beta(n,  k+1);
            b = log(Pi/4) + beta(n+2,k+1);
            
            if mask(k) == 1    
               c = log(Pt) + log((1-2*Ps)*(curr_marker_bit == y(n)) + Ps) + beta(n+1,k+1); 
            else
               c = log(Pt) + log(curr_Px_1*((1 == y(n))*(1-2*Ps) + Ps) + ...
                    (1-curr_Px_1)*((0 == y(n))*(1-2*Ps) + Ps)) + beta(n+1,k+1);  
            end
            beta(n,k) = max_star(a,b,c);
        

        elseif n == R
            % for state n = R 
            a = log(Pd) + beta(R,k+1);
            if mask(k) == 1  
                b = log(Pt) + log((1-2*Ps)*(curr_marker_bit == y(R)) + Ps) + beta(R+1,k+1);
            else
                b = log(Pt) +  log(curr_Px_1*((1 == y(n))*(1-2*Ps) + Ps) + ...
                    (1-curr_Px_1)*((0 == y(n))*(1-2*Ps) + Ps)) + beta(R+1,k+1); 
            end
            beta(R,k) = max_star(a,b);
        
        elseif n == R+1
            % for state n = R+1
            beta(R+1,k) = log(Pd)+beta(R+1,k+1);
        end
    end
 
end

end