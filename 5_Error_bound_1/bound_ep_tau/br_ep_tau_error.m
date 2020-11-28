function [err_bi, err_ep_tau_bound, err_Bhat] = ...
    br_ep_tau_error(B, A, N_hi, n, psi_ref, psi_low, c_low, sigma, r, p, xi_low, pc_solver)
% Compute bi-fidelity error and error bound quantities

% Inputs: 
% B - low-fidelity samples
% A - highfidelity samples
% N_hi - number of high-fidelity samples
% n - number of samples used to compute bound
% psi_ref - pc basis of complete high fidelity samples
% psi_low - pc basis of complete low fidelity samples
% sigma - tolerance on residual in spgl1 solver
% c_low - low fidelity pc coefficients
% r - truncation of KL expansion 

% Outputs: 
% err_bi - bi-fidelity error
% err_ep_tau_bound - bi-fidelity epsilon tau error bound

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample = datasample(1:size(A,1), N_hi, 'Replace', false); 
sample = 1:N_hi; 

%%% High fidelity model - limited samples N_hi
% xi_hi = xi_ref(sample,:); 
u_hi = A(:,sample)'; 
psi_hi= psi_ref(sample,:); 

opts = spgSetParms('iterations',10000,'verbosity',0);

%             Hi fid PCE solved for via ell_1 minimization
c_hi = spg_mmv(psi_hi,u_hi,sigma*norm(A),opts);            
%             c_hi = psi_hi\u_hi; 

c_hi = c_hi';
            
[~, ~, c_bi,~, ~, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma); 
            
%%% Bi-fidelity error calculation
psi_bi = psi_ref(:,2:end)*alpha2';
%psi_bi_ref = psi_ref(:,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi = [ones(size(A, 2), 1) psi_bi]; 
% remove r+1 column before solving for coefficients.
psi_bi = psi_bi(:,1:end); 
A_hat = (psi_bi*c_bi')'; 
err_bi = norm(A_hat-A)/norm(A); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(A,2); 

%%% Step 1: epsilon(tau)
% similar to mat_id_est_one_normal.m

% % Obtain column skeleton of P
% [P_s,ix] = matrixIDvR(B,r);

% Error bound inputs
% normC = norm(P_s);
sb = svd(B); 

%%% Bi-fidelity - With B - transposed 
%%% Figure out how to include this elegantly - ie compute early in advance,
%%% only hand error to the bound function. 

[c_low_L, ~] = my_pce(xi_low, p, B', sigma, pc_solver, 0); 
[~, ~, c_bi_L,...
                ~, ~, alpha_L]= ... 
                BR_FN(B',psi_low,c_low_L,c_low_L,r,sigma); 

psi_bi_L = psi_low(:,2:end)*alpha_L';
%psi_bi_ref = psi_ref(:,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi_L = [ones(size(B, 2), 1) psi_bi_L]; 
B_hat = (psi_bi_L*c_bi_L')'; 
err_Bhat = norm(B_hat-B); 

% Subset of vectors for bi-fidelity error estimate 
% Ensure sample for BR used + additional samples to reach
% total of n
rand_sample = [sample, getfield(setxor(sample,1:N), {1:n-numel(sample)})];

B_R = B(:,rand_sample);
A_R = A(:,rand_sample);

[err_ep_tau_bound, ~] = ...
    br_ep_tau_bound(B_R, A_R, err_Bhat, sb, n, N);

end

