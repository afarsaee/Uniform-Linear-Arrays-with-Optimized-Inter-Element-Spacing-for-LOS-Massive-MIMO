clc
clear all
close all
%% p2 over inter-element spacing d
addpath('func');
flag_write = 1;
M_antenna = 64;
power_vec = 0:0.0025:2.33; % to have range of spacing showed in Fig. 4
d_spacing = 0.5 * (2.^power_vec);  
half_FoV_deg = 90;
psi_thr = (1/(M_antenna)) * (0.5./d_spacing);
T = 1./d_spacing;
%% Main loop
n_sim = length(d_spacing);
n_sample = n_sim;
p2_spacing = zeros(1,n_sim);
for i_sim = 1:n_sim
   p2_spacing(i_sim) = find_p2_new(half_FoV_deg, T(i_sim), psi_thr(i_sim));
end
%% plotting the results
figure;
plot(d_spacing,p2_spacing)
title('p2 over spacing,     for \rho_{thr} = 2/\pi');
xlabel('d spacing');
ylabel('probability');
%% Writing the results
if flag_write == 1
    name = sprintf('p2_over_d_spacing_%d_%d.txt',2*half_FoV_deg,M_antenna);
    fp2       = fopen(name,'w');
    n_write = length(p2_spacing);
    step_size = floor(n_write/n_sample);
    for i = 1:step_size:n_write
       fprintf(fp2,'%0.6f %2.6f\n', d_spacing(i) ,p2_spacing(i));
    end
    fclose(fp2);
end    