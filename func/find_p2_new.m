 function p2_final = find_p2_new(half_FoV_deg, T, psi_thr)
 %% initialization
phi_min_deg = 90 - half_FoV_deg;
phi_max_deg = 90 + half_FoV_deg;

FoV2 = (pi*((2*half_FoV_deg)/180))^2;
phi_min_rad = (phi_min_deg/180)*pi;
phi_max_rad = (phi_max_deg/180)*pi;
%% part 1, alpha0
funA = @(x) acos(cos(x)-psi_thr) - acos(cos(x)+psi_thr);
qA = real(integral(funA,acos(cos(phi_min_rad) - psi_thr),acos(cos(phi_max_rad) + psi_thr)));
p2_main = qA/FoV2;
%% part 2, alpha0
funA = @(x) acos(cos(x)-psi_thr) - phi_min_rad;
qA_part2 = real(integral(funA,phi_min_rad,acos(cos(phi_min_rad) - psi_thr)));
p2_side = (qA_part2)/FoV2;
%% part 3, alpha0
funA = @(x) phi_max_rad - acos(cos(x)+psi_thr);
qA_part3 = real(integral(funA,acos(cos(phi_max_rad) + psi_thr),phi_max_rad));
p2_side_3 = (qA_part3)/FoV2;
%% alpha_infty
% #period that are located in the range of psi
% n_period = ceil((2/T))-1;

% Gamma_min = min((n_period+1)*T-psi_thr,2*cosd(phi_min_deg));
% Gamma_max = min((n_period+1)*T+psi_thr,2*cosd(phi_min_deg));
% p2_side_infty = 0;
% if Gamma_max - Gamma_min >= 0
%     funA = @(x) acos(cos(x)+Gamma_min) - acos(cos(x)+Gamma_max);
%     
%     qA_part3 = real(integral(funA,acos(cos(phi_min_rad) -Gamma_min),phi_max_rad));
%     p2_side_infty = (2*qA_part3)/FoV2;
% end
%% alpha_i
p2_main_T = 0;
i = 1;
flag_B = 0;
while i*T - psi_thr < 2*cos(phi_min_rad)
    % if the whole peak locates insider the psi interval
    if i*T + psi_thr < 2*cos(phi_min_rad)
        funA = @(x) acos(cos(x)-psi_thr + i*T) - acos(cos(x)+psi_thr + i*T);
        arg_max = phi_max_rad;
        arg_min = acos(cos(phi_min_rad) - psi_thr - i*T);
        
        flag_B = 1;
        funB = @(x) acos(cos(x)-psi_thr + i*T) - phi_min_rad;
        arg_maxB = acos(cos(phi_min_rad) - i*T - psi_thr);
        arg_minB = acos(cos(phi_min_rad) - i*T + psi_thr);
    elseif i*T + psi_thr >= 2*cos(phi_min_rad)
        % for the alpha_infty:
        funA = @(x) acos(cos(x)-psi_thr + i*T) - phi_min_rad;
        arg_max = phi_max_rad;
        arg_min = acos(cos(phi_min_rad) + psi_thr - i*T);
        flag_B = 0;
    end
    
    if arg_max > arg_min
        qA = integral(funA,arg_min,arg_max);
        if imag(qA) ~= 0
            error('rrl');
        end
    else
        qA = 0;
        error('not possible')
    end
    
    if flag_B == 1
        if arg_maxB > arg_minB
            qB = (integral(funB,arg_minB,arg_maxB));
            if imag(qB) ~= 0
                error('rrl');
            end
            qA = qA + qB;
        else
            error('not possible')
        end
    end
    
    p2_main_T = p2_main_T + (qA/FoV2);
    i = i + 1;
end
%%
ptot = p2_side_3 + p2_side + p2_main;
ptot_T = 2*(p2_main_T);
p2_final = ptot + ptot_T;
end