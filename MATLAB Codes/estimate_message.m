function [m_est, llr] = estimate_message(y, alpha, beta, Pd, Pi, Ps, T)
    
    % Total number of received bits
    R = length(y);

    % Likelihood arrays for individual bits
    l0 = zeros(1, T); 
    l1 = zeros(1, T);
    
    Pt = 1 - Pi - Pd;
    
    for k = 1:T

        % Get alpha vector corresponding to alpha vector @ stage k
        alpha_k = alpha(:,k);
        % Get beta vector corresponding to beta vector @ stage (k+1)
        beta_k = beta(:,k+1);
        % init the ins, del, trans0 and trans1 vectors
        ins = -realmax;
        del = -realmax;
        trans0 = -realmax;
        trans1 = -realmax;

        for n = max(1,R+1-2*(T-k+1)) : min(2*(k-1)+1, R+1)
            
            % vectors related to others
            if n == R + 1
               
                del = max_star(del, log(Pd)+beta_k(n)+alpha_k(n));
                
            elseif n==R
                
                del = max_star(del, log(Pd)+beta_k(n)+alpha_k(n));
                trans0 = max_star(trans0, log(Pt)+beta_k(n+1)+alpha_k(n)+log((1-Ps)*(y(n)==0) + Ps));
                trans1 = max_star(trans1, log(Pt)+beta_k(n+1)+alpha_k(n)+log((1-Ps)*(y(n)==1) + Ps));           
            else 
                ins = max_star(ins, log(Pi/4)+beta_k(n+2)+alpha_k(n));
                del = max_star(del, log(Pd)+beta_k(n)+alpha_k(n));
                trans0 = max_star(trans0, log(Pt)+beta_k(n+1)+alpha_k(n)+log((1-Ps)*(y(n)==0) + Ps));
                trans1 = max_star(trans1, log(Pt)+beta_k(n+1)+alpha_k(n)+log((1-Ps)*(y(n)==1) + Ps));
            end
            
          
        end
        l0(k) = max_star(trans0, del, ins);
        l1(k) = max_star(trans1, del, ins);
    end
   
   % Calculate LLR
   llr = l0 - l1;
  
 
   % message estimates
   m_est = llr < 0;
   

end