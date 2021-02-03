clc
clear all
close all
%%
tic
flag_write = 1;     %writing the results
rng('default');
addpath('func\');
%% Simulation Config.
n_realization = 100000;     % # of runs
n_user = 2:32;         % # of users
M_antenna = 64;   % # of antennas at the BS
% change the following spacing to get curves in Fig. 5
d_alpha = 0.5;      % the normalized spacing (d/lambda) --> 0.5, 0.994815 2.492
T = 1/d_alpha;      % T is used for finding the threshold
psi_thr = (1/(M_antenna))*(T/2);
rho_psi_thr = (sin(M_antenna*pi*(psi_thr*d_alpha)))/(M_antenna*sin(pi*(psi_thr*d_alpha))); % M*(2pi/lambda)*(d/2)*psi
two_over_pi = rho_psi_thr; % |\rho| > 2/pi --> highly correlated
phi_min =  0;    % min angle for the field-of-view (FoV)
phi_max =  180;   % max angle for the FoV
half_FoV_deg = (phi_max-phi_min)/2;
%% Main Loop for different number of users
n_sim = length(n_user);
p_at_least_two_cor_measured = zeros(1,n_sim);
for i_sim = 1:n_sim
    phi_ref = phi_min + (phi_max-phi_min).*rand(n_user(i_sim),n_realization);
    check_prob_all = 0;
    %% Main loop for n_user(i_sim)
	for i_realization = 1:n_realization
        %% First find rho_ij
        % read the phi values
        phi = phi_ref(:,i_realization);
        % find psi values, put 1 for the diagonal elements
        psi = ones(n_user(i_sim));
        for i = 1:n_user(i_sim)
            for j = 1:n_user(i_sim)
                if j ~= i
                    psi(i,j) = cosd(phi(i)) - cosd(phi(j));
                end
            end
        end
        % find the rho function sin(M*pi*d/lambda*psi)/sin(pi*d/lambda*psi)
        numerator_rho_func_sin   =   M_antenna*pi*d_alpha*psi;
        denumerator_rho_func_sin =   pi*d_alpha*psi;

        numerator_func_rho   = sin(numerator_rho_func_sin);
        denumerator_func_rho = sin(denumerator_rho_func_sin);
        % rho = 1/M (...)
        rho_ij_values = abs((1/M_antenna)*(numerator_func_rho./denumerator_func_rho));
        rho_ij_values = rho_ij_values - eye(n_user(i_sim));
        % check wether there is a problem (zero division)
        if sum(sum(isnan(rho_ij_values))) ~= 0
            error('error');
        end
        lower_part = tril(rho_ij_values);
        rho_lower_part = lower_part;
        %% look at the correlation and find any users that has 
        % a rho higher than 2/pi
        len_temp = length(find(rho_lower_part > two_over_pi));
        [a,b] = find(rho_lower_part > two_over_pi);
        array_ab = [a,b];
        len_distinguished_users = length(unique(array_ab));
        if len_temp >= 1
           check_prob_all = check_prob_all + 1;
        end
    end
    p_at_least_two_cor_measured(i_sim) = check_prob_all/n_realization;
end
%% p_k estimation and measured
figure;
plot(n_user,p_at_least_two_cor_measured);
title(sprintf('probability having at least two correlated users, #antennas = %d',M_antenna));
xlabel('number of users');
ylabel('probability');
%% Writing the results
if flag_write == 1
    fname_at_least_two_measured = sprintf('prob_correlated_at_least_two_measured_%d_another_figure.txt',M_antenna);
	fprob_at_least_two_measured = fopen(fname_at_least_two_measured,'w');
    
    n_write = length(p_at_least_two_cor_measured);
    for i = 1:n_write
       fprintf(fprob_at_least_two_measured,'%0.6f %2.8f\n', n_user(i) , p_at_least_two_cor_measured(i));
    end
    fclose(fprob_at_least_two_measured);
end
%%
toc