function [sum_rate_out, SNR_out] = find_ZF_SINR_max_min(H_in, Ptot)
    [n_user, ~] = size(H_in);
    UZF = normalize_columns(pinv(H_in));
    Gammai = abs(diag(H_in*UZF));
    one_over_Gammai = 1./Gammai;
    sum_denom_ZF_all = sum(one_over_Gammai.^2);
    SNR_out = Ptot/sum_denom_ZF_all;
    sum_rate_out = n_user * log2(1+SNR_out);
end