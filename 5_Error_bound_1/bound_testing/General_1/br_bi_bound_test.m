function [err_bi, err_bound, ahat_error_est, err_bound_N, err_bound_n,...
    err_Bhat, err_T, err_E, omega_term, q_N, q_n, Ir, ...
    err_E_exact,err_T_exact,err_CS_1,err_CS_2, A_hat, err_bi_vec_mean, err_bi_vec_med] = ...
    br_bi_bound_test(B, A, n, R, psi_ref, c_low, sigma, r)
% Compute bi-fidelity error and error bound quantities

% Inputs: 
% B - low-fidelity samples
% A - highfidelity samples
% n - number of high-fidelity samples to find bi basis
% R - number of high-fidelity samples to compute bound
% psi_ref - pc basis of complete high fidelity samples
% psi_low - pc basis of complete low fidelity samples
% sigma - tolerance on residual in spgl1 solver
% c_low - low fidelity pc coefficients
% r - truncation of KL expansion 

% Outputs: 
% err_bi - bi-fidelity error
% err_bound - bi-fidelity error bound
% ahat_error_est - Deterministic bound component
% err_bound_N - N samples component
% err_bound_n - n samples component
% err_Bhat - error estimating low-fidelity data - B 
% err_T - error component with T estimated
% err_E - error component with E estimated
% omega_term - max value of row in B term
% q_N - q for N samples
% q_n - q for n samples
% Ir - index of rank that optimizes epsilon tau bound
% err_T_exact - error component with T computed exactly
% err_E_exact - error component with E computed exactly
% err_CS_1 - E error if bound calculated pre-cauchy schwartz
% err_CS_2 - T error if bound calculated post-cauchy schwartz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(A,2); 

n_sample = datasample(1:N, n, 'Replace', false); 
% sample = 1:n; 

%%% High fidelity model - limited samples n
% xi_hi = xi_ref(sample,:); 
u_hi = A(:,n_sample)'; 
psi_hi= psi_ref(n_sample,:); 

opts = spgSetParms('iterations',10000,'verbosity',0);

%             Hi fid PCE solved for via ell_1 minimization
c_hi = spg_mmv(psi_hi,u_hi,sigma*norm(u_hi),opts);            
%             c_hi = psi_hi\u_hi; 

c_hi = c_hi';
            
[~, ~, c_bi,~, psi_bi_n, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma); 
            
%%% Bi-fidelity error calculation - for all N data samples
psi_bi_N = psi_ref(:,2:end)*alpha2';
%psi_bi_ref = psi_ref(:,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi_N = [ones(size(A, 2), 1) psi_bi_N]; 
% remove r+1 column before solving for coefficients.
psi_bi_N = psi_bi_N(:,1:end-1); 
A_hat = c_bi*psi_bi_N'; 
err_bi = norm(A_hat-A);
% err_bi_vec = vecnorm(A_hat-A,2, 2);
% e_check = A_hat-A; 
err_bi_vec_mean = mean(abs(A_hat-A),2);
err_bi_vec_med = median(abs(A_hat-A),2);
1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Step 1: epsilon(tau)
% similar to mat_id_est_one_normal.m

%%% Error bound inputs
sb = svd(B); 

% Bi-fidelity - With B - transposed 
% A_hat = c_bi*psi_bi'; 
B_hat = B*pinv(psi_bi_N')*psi_bi_N';
err_Bhat = norm(B_hat-B); 

% check norm: norm(pinv(psi_bi')*psi_bi')

% Is this Bhat correct? Using this bi-fidelity method where pc is limited
% approximation. 

% Subset of vectors for bi-fidelity error estimate 
% Ensure sample for BR used + additional samples to reach
% total of R if necessary. 

if R >= n 
    R_sample = [n_sample, datasample(setxor(n_sample,1:N), ...
        R-numel(n_sample), 'Replace', false)]; 
else
    R_sample = datasample(n_sample, R, 'Replace', false); 
end

B_R = B(:,R_sample);
A_R = A(:,R_sample);

A_n = u_hi';

[err_bound, ahat_error_est, err_bound_N, err_bound_n, ... 
err_T, err_E, omega_term, q_N, q_n, Ir] = ...
    br_bound(B_R, A_R, A_n, err_Bhat, sb, R, N, n, psi_bi_n);

% err_bound = 2*err_E +err_T*err_Bhat + 8*N^0.5*omega_term*(N^(-q_N)+n^(-q_n));

[err_E_exact,err_T_exact,err_CS_1,err_CS_2] = br_E_T(B,A, Ir, psi_bi_N');

% err_bound_ET = 2*err_E+err_T*err_Bhat; 
% E_check_rel = (err_E - err_E_exact)/err_E;
% T_check_rel = (err_T - err_T_exact)/err_T;
% CS_check_rel = (err_CS_2-err_CS_1)/err_CS_2;

end

