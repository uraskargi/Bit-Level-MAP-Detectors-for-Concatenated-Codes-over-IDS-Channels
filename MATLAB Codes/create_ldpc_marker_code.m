function [c, mask] = create_ldpc_marker_code(m, Nc, marker_code, perm, start_with_marker)
    
    % Params - marker code
    Nr = length(marker_code); 
    N = Nr + Nc; 
    rm = Nc/(Nc+Nr); 
    
    % ldpc encode - we are using only zero codeword
    % If Nc does not divide codeword size, then add zero padding
    if  rem(length(m), Nc)~= 0
        x = [m randi([0,1],1,Nc-rem(length(m), Nc))];
    else
        x = m;
    end

    % Permute the bits
    x = x(perm);
    xtotal = length(x);
    
    if start_with_marker == 0
        % create marker coded bit
        c = zeros(1, xtotal/rm); 
        % marker bit
        mask = zeros(1, xtotal/rm);

        for i = 1:xtotal/Nc
            low_ind = N*(i-1)+1; % low ind
            high_ind = N*i; % high ind
            low_ind_m = (N-Nr)*(i-1) + 1;
            high_ind_m = (N-Nr)*i;
            c(1, low_ind : high_ind - Nr) = x(1, low_ind_m: high_ind_m);
            c(1, high_ind - Nr + 1: high_ind) = marker_code;
            mask(1, high_ind - Nr + 1: high_ind) = ones(1, Nr);
        end
    else
        % create marker coded bit
        c = zeros(1, xtotal/rm + Nr); 
        c(1,1:Nr) = marker_code;
        % marker bit
        mask = zeros(1, xtotal/rm + Nr);
        mask(1,1:Nr) = ones(1, Nr);

        for i = 1:xtotal/Nc
            low_ind = N*(i-1)+1; % low ind
            high_ind = N*i; % high ind
            low_ind_m = (N-Nr)*(i-1) + 1;
            high_ind_m = (N-Nr)*i;
            c(1, Nr + low_ind : high_ind) = x(1, low_ind_m: high_ind_m);
            c(1, high_ind  + 1: high_ind + Nr) = marker_code;
            mask(1, high_ind + 1: high_ind + Nr) = ones(1, Nr);
        end
    end

    

    
end