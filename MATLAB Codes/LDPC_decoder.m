function [m_est,llr] = LDPC_decoder(c_llr, H, iter_num, early_finish)
%
%   c_llr = channel llr values
%   H = parity check matrix 
%   iter_num = total iteration number
%

size_H = size(H);
n = size_H(2);
m = size_H(1);
% initiliaze the messages
VN_to_CN = H.*c_llr;
CN_to_VN = zeros(size_H);

for iter = 1:iter_num
    CN_to_VN = zeros(size_H);

    % Calculate VN to CN messages
    for i = 1:m
        
        % Get current row vector of H
        curr_row = H(i,:);
        % Find where this equals to non-zero
        ind = find(curr_row);
        % get messages coming to CN
        mes = VN_to_CN(i,ind);
        
        for j = 1:length(ind)
    
            % Delete the corresponding CN_to_VN meesage where we send info
            mes_upd = mes;
            mes_upd(j) = [];
            
            % Find current CN to VN message
            CN_to_VN(i, ind(j)) = 2*atanh(min(0.9999999999999, ...
                max(-0.9999999999999, prod(tanh(1/2*mes_upd)))));
        end
    end
    
    % Estimate the messages
    llr = sum(CN_to_VN, 1) + c_llr;
    m_est = llr < 0;
    %ext_llr = llr - c_llr;
    
    if early_finish == 1
        if sum(mod(m_est * H.',2), 'all') == 0
            %disp('Decoding is finished at itetation')
            %disp(iter)
            break;
        end
    end

    VN_to_CN = zeros(size_H);
    % Calculate VN to CN messages
    for j = 1:n
        
        % Get current column vector of H
        curr_col = H(:,j);
        % Find where this equals to non-zero
        ind = find(curr_col).';
        % get messages coming to CN
        mes = CN_to_VN(ind, j);
        
        for i = 1:length(ind)
    
            % Delete the corresponding CN_to_VN meesage where we send info
            mes_upd = mes;
            mes_upd(i) = [];
            
            % Find current VN to CN message
            VN_to_CN(ind(i),j) = sum(mes_upd)+c_llr(j);
        end
    end
end

% Estimate the messages
llr = sum(CN_to_VN, 1) + c_llr;
m_est = llr <= 0;
%llr = llr - c_llr;


end





