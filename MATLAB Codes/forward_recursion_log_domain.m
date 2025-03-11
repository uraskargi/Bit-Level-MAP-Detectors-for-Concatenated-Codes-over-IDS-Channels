function alpha = forward_recursion_log_domain(y,llr_apriori,T,Pd,Pi,Ps,mask,marker_code,N,Nc)
%FORWARD_RECURSION_LOG_DOMAIN Calculates forward probabilities for bit-level
%   marker code BCJR coding in the log domain.
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
%       alpha       - Forward probabilities in the log domain
%
%   This function implements the forward recursion in the log domain for
%   bit-level marker code BCJR coding. It calculates the forward
%   probabilities (alpha) for each received bit position and transmitted
%   bit position, taking into account the a priori information, deletion,
%   insertion, and substitution events.

% Get a-priori probs from LLR information
Px_1 = 1./(1+exp(llr_apriori));     % P(x_i = 1) for i = 1,2,...,n

R = length(y);                      % Recieved bit number
alpha = -realmax*ones(R+1,T+1);     % Forward recursion inits
alpha(1,1) = 0;                     % Init. the starting point
ind = 1;                            % Denotes the index of the message bit
ind_marker = 1;                     % Denotes the index of the marker
Pt = 1-Pd-Pi;                       % Transmission probability
Nr = N- Nc;

% From stage 1 to T
for k = 1:T 
    
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
        curr_marker_bit = marker_code(ind_marker);  
        ind_marker = ind_marker + 1;
    end

    for n = max(1,R+1-2*(T-k)):min(2*k+1,R+1)
        
        % For state n = 1, only deletion event occurs
        if n == 1 
            alpha(n,k+1) = log(Pd) + alpha(n,k);  
        % For state n = 2, deletion or transmission events can occur
        elseif n == 2
            a = log(Pd) + alpha(2,k);
            % If this is market bit position, we need the marker bit @ that position
            if mask(k) == 1      
                b = log(Pt) + log((1-2*Ps)*(curr_marker_bit == y(n-1)) + Ps) + alpha(n-1,k);
            else
                b = log(Pt) + log(curr_Px_1*((1 == y(n-1))*(1-2*Ps) + Ps) + ...
                    (1-curr_Px_1)*((0 == y(n-1))*(1-2*Ps) + Ps)) + alpha(n-1,k);
            end
            alpha(n,k+1) = max_star(a,b);
        
        % For states bigger, deletion, transmission
        % or insertion can occur
        elseif n > 2 

            a = log(Pd)   + alpha(n,k);
            b = log(Pi/4) + alpha(n-2,k);
            if mask(k) == 1
               c = log(Pt) + log((1-2*Ps)*(curr_marker_bit == y(n-1)) + Ps) + alpha(n-1,k); 
            else
               c = log(Pt) + log(curr_Px_1*((1 == y(n-1))*(1-2*Ps) + Ps) + ...
                    (1-curr_Px_1)*((0 == y(n-1))*(1-2*Ps) + Ps)) + alpha(n-1,k);
            end 
            alpha(n, k+1) = max_star(a,b,c);
        end

    end

end
end




