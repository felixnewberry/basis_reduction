function [err_bi, err_ep_tau_bound, err_Bhat, err_E, err_T, Ir, err_E_exact, ...
    err_T_exact, err_CS_pre, err_CS_post, err_bi_N, err_bi_inf_N, err_bi_inf_n] = ...
    br_ep_tau_error_dev_Nsubset(B, A, N_hi, n, psi_ref, psi_low, c_low, sigma, r, p, xi_low, pc_solver, n_est)
% Compute bi-fidelity error and error bound quantities

% Inputs: 
% B - low-fidelity samples
% A - highfidelity samples - first n_est correspond to B - then additional
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

sample = datasample(1:n_est, N_hi, 'Replace', false); 
% sample = 1:N_hi; 

%%% High fidelity model - limited samples N_hi
% xi_hi = xi_ref(sample,:); 
u_hi = A(:,sample)'; 
psi_hi= psi_ref(sample,:); 

% opts = spgSetParms('iterations',10000,'verbosity',0);

% pc_val = 0; 
% pc_solver = 2; 
% [c_hi, psi_ref, ~] = my_pce(xi_ref, p, u_hi, sigma, pc_solver, pc_val); 

%  = spg_mmv(psi_hi,u_hi,sigma*norm(u_hi),opts);            
% %             c_hi = psi_hi\u_hi; 
% c_hi = c_hi';

weights = get_matrix_weights(psi_hi);
Psiweights = psi_hi*weights;    
norm_u_vec = vecnorm(u_hi,2); 
%                     c_data{n_points} = []; 
n_points = size(u_hi,2); 
P = size(psi_hi,2); 
% c_hi = zeros(n_points, P);
c_data{n_points} = []; 
parfor i_points=1:n_points
    opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

    delta = sigma*norm_u_vec(i_points);
    c_data{i_points}.c_vec = weights*spg_bpdn(Psiweights,u_hi(:,i_points),delta,opts);
end

c_hi = zeros(n_points, P);
for i_points=1:n_points
    c_hi(i_points,:) = c_data{i_points}.c_vec; 
end
%%% High fidelity coeficients are used only to compute eigen values. Don't
%%% do spg in interests of time. 
% n_points = size(A,1); 	
% P = size(psi_hi,2); 
% 
% weights = get_matrix_weights(psi_hi);
% Psiweights = psi_hi*weights;
% 
% norm_u_vec = vecnorm(u_hi,2); 
% 
% c_hi = zeros(n_points, P);
% for i_points=1:n_points
%     opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);
% 
%     delta = sigma*norm_u_vec(i_points);
%     c_hi(i_points,:)  = weights*spg_bpdn(Psiweights,u_hi(:,i_points),delta,opts);
% end
                    
            
[~, ~, c_bi,~, ~, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma); 
            
%%% Bi-fidelity error calculation
psi_bi = psi_ref(1:n_est,2:end)*alpha2';
%psi_bi_ref = psi_ref(:,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi = [ones(size(A(:,1:n_est), 2), 1) psi_bi]; 
% remove r+1 column before solving for coefficients.
psi_bi = psi_bi(:,1:end); 
A_hat_n = (psi_bi*c_bi')'; 
err_bi = norm(A_hat_n-A(:,1:n_est));


%%% estimate with all samples

[~, ~, c_bi_N,~, ~, ~]= BR_FN(A(:,1:n_est)',psi_ref(1:n_est,:),c_hi,c_low,r,sigma); 
A_hat_N = (psi_bi*c_bi_N')'; 
err_bi_N = norm(A_hat_n-A_hat_N);

[~, ~, c_bi_inf,~, ~, ~]= BR_FN(A',psi_ref,c_hi,c_low,r,sigma); 
A_hat_inf = (psi_bi*c_bi_inf')'; 
err_bi_inf_n = norm(A_hat_n-A_hat_inf);
err_bi_inf_N = norm(A_hat_N-A_hat_inf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = n_est; 

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

[err_ep_tau_bound, err_E, err_T, Ir] = ...
    br_ep_tau_bound_dev(B_R, A_R, err_Bhat, sb, n, N);

[err_E_exact,err_T_exact, err_CS_pre, err_CS_post] = br_E_T_dev(B, A(:,1:n_est), Ir, psi_bi_L');

end

