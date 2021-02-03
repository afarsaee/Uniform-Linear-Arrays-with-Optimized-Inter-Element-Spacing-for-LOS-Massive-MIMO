%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functions outputs an uplink channel %
% for LOS, it supports ULA and UPA (square array)
% Users are distributed uniformly in R^2, phi, and theta :)
                        %% !!! Important !!!
                % I'm not sure how we should model the NLOS !
                % I mean, I can take the pathloss from 3GPP file, 
                % but what is NLOS channel? is it Rayleigh?
%% Inputs:
% M:                                #antennas at the BS
% K:                                #users (single-antenna)
% f0_gig:                           carrier frequency in hz
% R_max:                            maximum radius of user
% R_min:                            minimum radius of user
% alpha_LOS:                        1 if LOS, 0 if Rayleigh
% theta_min:                        minimum phi_angle (horizontal)
% theta_max:                        maximum phi_angle (horizontal)
% min_spacing_distance_user:        minimum allowable distance between users
% spacing_array:                    0.5 if half-wavelength spacing
% shadowing:                        1 if you want shadowing (based on 3GPP)
% flag_2D:                          1 if you want a uniform planar array (square)
%% Outputs:
% channel_unit_norm:                uplink channel with unit norm column vectors
% channel_gain:                     channel norms of for the users
% path_loss_dB:                     path loss (channel norms squared) for each TX-RX antennas
% channel_nonnormalized:            non-normalized uplink channel matrix

function [H_out] = gen_LOS_channel(M_ant,K_user,phi_min, phi_max, min_spacing_phi, spacing_array)
%% Read the inputs
n_user = K_user;             % # single-antenna user
n_bs = M_ant;               % # bs antennas
H_out       = zeros(n_bs,n_user);
%% if ULA:
%% Distributing the users uniformly in a sector
% distributed uniformly in R^2, and uniformly in phi
% and uniformly in theta
flag_user_spacing = 1;  % flag is true as long as the users do not meet the min spacing
while flag_user_spacing == 1
    % phi uniformly distributed \in [\phi_min,\phi_max]
    phi_user_deg = phi_min + (phi_max-phi_min).*rand(1,n_user);
    phi_user_deg_sorted = sort(phi_user_deg);
    for i_user = 1:n_user
        for j_ant = 1:n_bs
            H_out(j_ant,i_user) = exp(1j * (j_ant-1) * ((2*pi)) * spacing_array * cosd(phi_user_deg(i_user))); 
        end
        H_out(:,i_user) = H_out(:,i_user)/sqrt(n_bs);
    end
    %% Check whether two users are very close to each other:
    % when both the dimensions of the users are closer than the min
    % spacing, then the users are really really close!
    % we avoid these scenarios! by repeating the generating the users
    phi_user_difference_flag =  abs(diff(phi_user_deg_sorted)) <= min_spacing_phi;
    if sum(phi_user_difference_flag) == 0
        % if no co-located users:
    	flag_user_spacing = 0;
    else
        % if at least two co-located users
    	flag_user_spacing = 1;
    end
end
end