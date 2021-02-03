clc
clear all
close all
tic
%% Figure 6, 
% CDF plot for ZF sum-rate, with no dropping
% using the designed ULAs compared to half-wavelength array
%%
rng('default');
addpath('func');         % adding the path for func
flag_write_CDF = 1;      % flag to write and run for CDF (one SNR)
%% LOS Configuration
% 100K channel, 64-antennas ULA serving 6 users
% phi \in 0:180
alpha_LOS = 1;
min_spacing_phi_deg = 0.01;
n_channel = 100000;  n_bs = 64;  n_user_ref = 6;
phi_min = 0; phi_max = 180; 
% Change the following spacing to get the curves in Fig. 6
spacing_array = 2.492;% 0.5, 0.994815 2.492
%% Simulation parameters
% /in Favorable propagation a data rate of 6 bit/transmission is achieved per user/
bits_orthogonal = 6;

% requried SNR
mySNRdB = 10*log10((2^bits_orthogonal)-1);
mySNR = 10.^(mySNRdB/10);
%% Defining the variables
sum_rate_ZF_full    = 0;
CDFSNR_ZF_full    = zeros(n_channel,1);
Ptot_ref_per_SNR = mySNR * n_user_ref;   
%% Repeat the simulation for n_channel realizations
for i_channel = 1:n_channel
   n_user_THP_Dropped = n_user_ref;
   [H_out_uplink] = gen_LOS_channel(n_bs,n_user_ref,phi_min, phi_max, min_spacing_phi_deg, spacing_array);
   channel_current_downlink = (H_out_uplink');
   Ptot = Ptot_ref_per_SNR;
   H_ZF = channel_current_downlink;
   %% ZF + Full
   [sum_rate_out_ZF_full, SINR_ZF_all] = find_ZF_SINR_max_min(H_ZF, Ptot);
   sum_rate_ZF_full = sum_rate_ZF_full + sum_rate_out_ZF_full;
   CDFSNR_ZF_full(i_channel) = SINR_ZF_all;
end
sum_rate_ZF_full      = sum_rate_ZF_full/n_channel;
%% Show the CDF plot
figure;
% ZF full
h4 = cdfplot(CDFSNR_ZF_full(:));
CDF_SNRw_ZF_full_val          = get(h4,'YData');
CDF_SNRw_ZF_full              = get(h4,'XData');
index_5_percentile = find(CDF_SNRw_ZF_full_val >= 0.05);
val_005_ZF_full_SNRw = CDF_SNRw_ZF_full(index_5_percentile(1));
legend('ZF Full');
display(['5 percentile ZF full      = ',num2str(val_005_ZF_full_SNRw)]);
title('SNR');
%% Show the CDF of sum-rate
figure;
CDFSumRate_ZF_full = n_user_ref * log2(1 + CDFSNR_ZF_full);
h44 = cdfplot(CDFSumRate_ZF_full(:));
CDF_SumRatew_ZF_full_val = get(h44,'YData');
CDF_SumRatew_ZF_full     = get(h44,'XData');
index_5_percentile = find(CDF_SumRatew_ZF_full_val >= 0.05);
val_005_ZF_full_SumRatew = CDF_SumRatew_ZF_full(index_5_percentile(1));
legend('ZF Full');
display(['5 percentile ZF full      = ',num2str(val_005_ZF_full_SumRatew)]);
title('Sum-rate');
%% writing the CDF plots
if flag_write_CDF == 1
    for i_dummy = 1:1
        name_ZF_full        = sprintf('CDF_ZF_Full_dB_%d_%d_%2.2d.txt',n_bs,n_user_ref, spacing_array);
        name_ZF_fullSumRate = sprintf('CDF_SumRate_ZF_Full_dB_%d_%d_%2.2d.txt',n_bs,n_user_ref, spacing_array);
                
        file_arg_ZF_full        = fopen(name_ZF_full,'w');       
        file_arg_ZF_fullSumRate = fopen(name_ZF_fullSumRate,'w');
        
        n_write = length(CDF_SumRatew_ZF_full);
        n_step  = floor(n_write/1000);
        for i = 3:n_step:n_write
           fprintf(file_arg_ZF_fullSumRate ,'%0.6f %2.6f\n', CDF_SumRatew_ZF_full(i) ,CDF_SumRatew_ZF_full_val(i));
        end
        fclose(file_arg_ZF_fullSumRate);        
        n_write = length(CDF_SNRw_ZF_full);
        n_step  = floor(n_write/1000);
        for i = 3:n_step:n_write
           fprintf(file_arg_ZF_full ,'%0.6f %2.6f\n', 10*log10(CDF_SNRw_ZF_full(i)) ,CDF_SNRw_ZF_full_val(i));
        end
        fclose(file_arg_ZF_full);
    end
end
a = toc