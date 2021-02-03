function [H_dropped, n_user_dropped] = Drop_user_ZF_fixed(H, n_drop)
   [n_user,M_ant] = size(H);
   HHH = find_rho_ij(H);
   n_user_new = n_user;
   %% drop n_drop users
   while n_user_new >= 2 && n_user_new > n_user - n_drop
        % find the correlated users indexes
        [~,index_cor] = max(HHH(:));
        [cor_user1,cor_user2] = ind2sub(size(HHH),index_cor);
        HHH_new = (abs(HHH));
        HHH_new(cor_user1,cor_user2) = 0;
        HHH_new(cor_user2,cor_user1) = 0;

        compare_cor_user_1 = max(HHH_new(cor_user1,:));
        compare_cor_user_2 = max(HHH_new(cor_user2,:));
        % the new channel matrix
        H_new = zeros(n_user_new-1,M_ant);
        % keep the user which is less correlated to the others
        if compare_cor_user_1 > compare_cor_user_2
            H_new(1,:)      = H(cor_user2,:);
        else
            H_new(1,:)      = H(cor_user1,:);
        end
        % rebuild the channel matrix
        index_new = 1;
        for i = 2:n_user_new-1
            while index_new == cor_user1 || index_new == cor_user2
                index_new = index_new + 1;
            end
            H_new(i,:) = H(index_new,:);
            index_new = index_new + 1;
        end
        % update the number of user and channel
        n_user_new = n_user_new - 1;
        H = H_new;
        HHH = find_rho_ij(H);
   end
   n_user_dropped = n_user_new;   
   H_dropped      = H;
end