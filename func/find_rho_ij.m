%% This function finds the normalized spatial correlation for the channel H
% Input:
    % H: input channel matrix,
% Output:
    % rho_ij: the normalized spatial correlation matrix 
    % (diagonal elements are made 0)
%%
function rho_ij = find_rho_ij(H)
[n_user,~] = size(H);
rho_ij = zeros(n_user);
for i = 1:n_user
    for j = 1:n_user
        if j~= i
            norm_Hi = norm(H(i,:));
            norm_Hj = norm(H(j,:));
            rho_ij(i,j) = abs((H(i,:)*H(j,:)')/(norm_Hi*norm_Hj));
        end
    end 
end
end