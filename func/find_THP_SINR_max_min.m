function [sum_rate_out, SNR_out] = find_THP_SINR_max_min(H_in, Ptot, flag_low)
    [n_user, ~] = size(H_in);
    if flag_low == 0
        H_in = VBLAST(H_in);   % reorder the users
    elseif flag_low == 1
        H_in = VBLAST_low(H_in);   % reorder the users
    end
    % Now do the QR decomposition
    HUplink = H_in';
    [~,R] = qr(HUplink);
    % find L
    L = R';
    % find 1/l(i,i)
    one_over_Lii_correlated_users = 1./diag(L);
    % sum_i 1/l(i,i)^2
    sum_one_over_lii2 = sum(one_over_Lii_correlated_users.^2);
    % beta_THP = Ptot/sum_i 1/l(i,i)^2 
    % THP with max-min power control SINR
    beta_THP2 = Ptot/sum_one_over_lii2;

    SNR_out = beta_THP2;

    sum_rate_out = n_user * log2(1+(SNR_out));
end