clc
clear all
close all
%% plotting rho over psi
n_points = 2000;
flag_write = 1;
M = 10;
c = 0.3;
f = 30;
lambda = c/f;
k = (2*pi)/lambda;
d1_lambda = 1*lambda;
d2_half = lambda/2;

phi_min = 0;
phi_max = 180;
%%
figure;
psi = linspace(-2,2,n_points);
func_rho = @(saw_in,d_in) abs((sin(M*pi*d_in*saw_in))./(sin(pi*d_in*saw_in)))/M;
rho_d_lambda      = func_rho(psi,d1_lambda/lambda);
rho_d_lambda_half = func_rho(psi,d2_half/lambda);
plot(psi,rho_d_lambda);
hold on;
plot(psi,rho_d_lambda_half);
legend('lambda','lambda half');
%%
if flag_write == 1
    name_lambda = sprintf('rho_%d_lambda_%d_%d.txt',d1_lambda/lambda,M,phi_max);
    name_half   = sprintf('rho_%1.1f_%d_%d.txt',d2_half/lambda,M,phi_max);
	fpsilambda       = fopen(name_lambda,'w');
    fpsilambdahalf   = fopen(name_half,'w');
    n_write = length(psi);
    for i = 1:n_write
       fprintf(fpsilambda,'%0.6f %2.6f\n',       psi(i) ,rho_d_lambda(i));
       fprintf(fpsilambdahalf,'%0.6f %2.6f\n',   psi(i) ,rho_d_lambda_half(i));
    end
    fclose(fpsilambda);
    fclose(fpsilambdahalf);
end