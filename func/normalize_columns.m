% =========================================================================
% Normalizing the columns of input matrix to have unit-norm column vectors
%   -- inputs:
%       - U: input complex matrix
%   -- outputs: 
%       - Uout: output complex matrix with unit-norm column vectors
% -------------------------------------------------------------------------
% (c) Oct. 2019 Ashkan Farsaei
% e-mail: a.farsaee@tue.nl
% =========================================================================
function Uout = normalize_columns(U)

% read the number of columns
[~,col_size] = size(U);

% initialize the output
Uout = zeros(size(U));

% normalize the columns of input matrix and save it at the output
for i = 1:col_size
    Uout(:,i) = U(:,i)/norm(U(:,i));
end

end