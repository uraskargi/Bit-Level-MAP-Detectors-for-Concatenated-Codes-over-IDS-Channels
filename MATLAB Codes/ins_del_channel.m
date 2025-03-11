function [y, trans] = ins_del_channel(c, Pd, Pi, Ps)
    
    len_c = length(c); % codeword length
    Pt = 1 - Pd - Pi; % prob of cor trans
    
    % get trans probability for each case
    trans = randsample(['c','d','i'], ...
        len_c, true, [Pt, Pd, Pi]);

    % recieved signal
    y = [];

    for i = 1: len_c
        
        % correct transmission
        if trans(i) == 'c'
            rand_m = rem(c(i) + randsample([0,1],1, true, [1-Ps, Ps]), 2);
            y = [y rand_m];

        elseif trans(i) == 'i'
            i_bits = randi([0 1],1,2);
            y = [y i_bits];
        end

    end

end
