function[mean_bi, var_bi, c_bi, lam_hi, psi_bi, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma)
% Performs basis rotation of a polynomial chaos expansion to determine it's
% compact representation
% 
% Background:
% 1) Approximate the vector-valued quantity of interest from the low-fid
% model in a PC basis using l_1 minimization. This is performed prior to
% calling BR_FN so low fidelity PCE coefficients are a required input. 
% 2) Write the low-fid solution in an exact KL expansion to get the
% description of the KL random variables in the PC used. Can be pictured as
% rotating the PC basis with the KL basis (both orthonormal basis) while
% removing many of the unncessary basis functions. 
% 3) Set up an l_1 regression problem for the high fidelity model
% approximation to aquire a new reduced basis c_bi. 
% 
% Inputs:
% u_hi - vector/matrix of evaluations of the high fidelity quantity of
% interest. 
% psi_hi - measurement matrix of the high fidelity quantity of interest
% c_hi - high fidelity PCE coefficients
% c_low - low fidelity PCE coefficients
% c_ref - reference PCE coefficients
% p - not used... could remove? 
% r - rank of the approximation (i.e. the number of basis terms you want to
% end up with from the bi-fidelity model 
% 
% Outputs:
% mean_bi - mean of bi-fidelity model
% var_bi - variance of bi-fidelity model
% c_bi - reduced basis
% lam_ref - vector of eigenvalues of the reference model
% lam_hi - vector of eigenvalues of the high fidelity model
% lam_low - vector of eigenvalues of the low fidelity model
% sigma - tolerance on residual used in ell_1 minimization
%
% Requirements: 
% spgl1 toolbox

%% Bi-fidelity 

% Covariance of the low-fid PC solution. Note that the index starts at 2.
% We exclude the first column of c_low here since the basis functions are 
% entirely ones corresponding to measurement matrix of column 1 means. 
% Add it in later 


% size c_low is n_grid points (9) by basis functions P
cov_low = c_low(:, 2:end) * c_low(:, 2:end)'; 

% Eigenvalue decomposition to get the eigenvalues and eigenvectors of the
% nodal covariance of the low-fidelity model

% compute r + 1 here, truncate after having identified sigma_k+1 for errro
% bound. 
[v_low, lam_low] = eigs(cov_low, r+1); 

% Projection of low-fid PC solution onto the eigenvectors of covariance
% yields size r (approximation rank) by P -1 
alpha2 = diag(1./sqrt(diag(lam_low)))*(v_low'*c_low(:, 2:end)); 

% New reduced basis (note: psi_hi is independent of low or high fid. Just
% depends on number of samples available. 
psi_bi = psi_hi(:,2:end)*alpha2';
%psi_bi_ref = psi_ref(:,2:end)*alpha2';

% New reduced basis including column of ones
psi_bi = [ones(size(u_hi, 1), 1) psi_bi]; 

% Calculate svd for error bound, can I use eigenvalue? 
% S = svd(psi_bi);
% sv_rp1 = S(r+1); 
% sv_rp0 = S(r); 

% remove r+1 column before solving for coefficients.
psi_bi = psi_bi(:,1:end-1); 

% Necessary paramters for spgl1 toolbox
opts = spgSetParms('iterations', 10000, 'verbosity', 0); 

% l_1 minimization to solve for new reduced basis coefficients

c_bi = spg_mmv(psi_bi, u_hi, sigma*norm(u_hi), opts); 
% c_bi = psi_bi\u_hi; 
c_bi = c_bi'; 

%%% Bi-fidelity
% Statistics
mean_bi = c_bi(:,1); 
var_bi = sum(c_bi(:,2:end).^2, 2); 

%%% Eigenvalue generation for the reference model and low-fid model
% Eigenvalues of the reference model - from this you can get an
% interesting, though not entirely necessary comparison

% % Nodal covariance of the reference model
% cov_ref = c_ref(:, 2:end)*c_ref(:,2:end)';

% % Eigenvalue decomposition of nodal covariance
% [~,lam_ref] = eigs(cov_ref,r);

% Nodal covariance of the high fidelity model
cov_hi = c_hi(:,2:end)*c_hi(:,2:end)'; 

% Eigenvalue decomposition of nodal covariance
[~, lam_hi] = eigs(cov_hi,r); 

% eigenvalues for hi, reference, and low fidelity
lam_hi = diag(lam_hi);
% lam_ref = diag(lam_ref);
% lam_low = diag(lam_low); 
% lam_low = lam_low(1:r); 

end



