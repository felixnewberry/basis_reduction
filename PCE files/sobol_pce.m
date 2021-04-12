function s = sobol_pce( pc_coeff, index_pc_j)

% This function computes PC-based Sobol index, i.e., s, of the jth random variable 
% index_pc_j is the jth column of index_pc

P = size(index_pc_j,1); 

% Compute the variance
var = sum(pc_coeff(2:P,1).^2,1);

k = 1;
for i = 1:P
    if index_pc_j(i,1) ~= 0
        s_index(k,1) = i;
        k = k + 1;
    end
end

s = 0;
for i=1:size(s_index,1)
    s = s + pc_coeff(s_index(i,1),1)^2;
end

s = s/var;

end

