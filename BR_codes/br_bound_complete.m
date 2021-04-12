function [err_bi_mean, err_bi_sum, mu, rho_k, zeta_i_1, zeta_i_2, ...
    zeta_N_1, MID_efficacy, eff_36, bound_34, ...
    p_33, p_35, R] = ...
    br_bound_complete(B, A, N_hi, n, psi_ref, c_low, sigma, r, n_est, t, C)
% Compute bi-fidelity error and error bound quantities

% Inputs: 
% B - low-fidelity samples - normalized by fro of first n_est. 
% A - highfidelity samples - first n_est correspond to B - then additional
% for some problems (i.e. Lid driven cavity has 2200)
% N_hi - number of high-fidelity samples
% n - number of samples used to compute bound
% psi_ref - pc basis of complete high fidelity samples
% sigma - tolerance on residual in spgl1 solver
% c_low - low fidelity pc coefficients
% r - truncation of KL expansion 
% n_est - how many samples are we estimating? i.e. size of A and A_hat
% t - bound via moments
% C - constant for bound via moments

% Outputs: 
% err_bi_sum    - bi-fidelity error squared sum
% err_bi_mean    - bi-fidelity error each point squared
% mu            - coherence
% rho_k         - bi-fid error rank revealing QR
% Y_Nh_vec      - Y for each point in M_H
% theta_vec     - theta for each point in M_H
% U_bar         - evaluation of norm in (27) - squared

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PC Bi-fidelity and High fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sample = datasample(1:n_est, N_hi, 'Replace', false); 
% sample = 1:N_hi; 

%%% High fidelity model - limited samples N_hi
u_hi = A(:,sample)'; 
psi_hi= psi_ref(sample,:); 

weights = get_matrix_weights(psi_hi);
Psiweights = psi_hi*weights;    
norm_u_vec = vecnorm(u_hi,2); 

n_points = size(u_hi,2); 
P = size(psi_hi,2); 

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

%%% Calculate PC reduced basis
[~, ~, c_bi,~, ~, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma); 
            
%%% Bi-fidelity error calculation
psi_bi = psi_ref(1:n_est,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi = [ones(size(A(:,1:n_est), 2), 1) psi_bi(:,1:end-1)]; 
A_hat = (psi_bi*c_bi')'; 

%%% calculate squared error for each sample and each point. 
err_bi_mat = (A_hat-A(:,1:n_est)).^2;
%%% LHS of inequalities (19), (20), (22) and (34)
err_bi_mean = mean(err_bi_mat,2); 
%%% LHS of inequalities (23), (25)-(27) and (36)
err_bi_sum = sum(err_bi_mean);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Coherence for N_hi samples
A1 = psi_bi(:,2:end); 
mu = max(sum(abs(A1(sample,:)).^2,2));

% Term 1 in (22), (23) and every other appearence of the bound. 
term_1 = (1+4*mu/N_hi);

%%% Estimate truncation error
% RHS of (11) for n=N r=r_max=size(c_low,1)
r_max = rank(c_low(:, 2:end) * c_low(:, 2:end)');
[~, ~, c_bi_ideal,~, ~, alpha2_ideal]= BR_FN(A',psi_ref,c_hi,c_low,r_max,sigma); 
psi_bi_ideal_R = psi_ref(:,2:end)*alpha2_ideal';
psi_bi_ideal_R = [ones(size(A(:,:), 2), 1) psi_bi_ideal_R(:,1:end-1)]; 
A_hat_ideal_R = (psi_bi_ideal_R*c_bi_ideal')'; 
% RHS of (11) for n=N and r=r
[~, ~, c_bi_ideal_r,~, ~, alpha2_ideal_r]= BR_FN(A',psi_ref,c_hi,c_low,r,sigma); 
psi_bi_ideal_r = psi_ref(:,2:end)*alpha2_ideal_r';
psi_bi_ideal_r = [ones(size(A(:,:), 2), 1) psi_bi_ideal_r(:,1:end-1)]; 
A_hat_ideal_r = (psi_bi_ideal_r*c_bi_ideal_r')'; 

err_N_R = (A_hat_ideal_R-A_hat_ideal_r);
delta_R_N_i2 = mean(err_N_R.^2,2);

%%% estimate truncation error 2. (first estimate is 2 orders of magnitude
%%% too large). 
% RHS of (11) for n=n and r=R
[~, ~, c_bi_n_R,~, ~, alpha2_n_R]= BR_FN(u_hi,psi_hi,c_hi,c_low,r_max,sigma); 
psi_bi_n_R = psi_ref(:,2:end)*alpha2_n_R';
psi_bi_n_R = [ones(size(A(:,1:n_est), 2), 1) psi_bi_n_R(1:n_est,1:end-1)]; 
A_hat_n_R = (psi_bi_n_R*c_bi_n_R')'; 
% RHS of (11) for n=n and r=r (already have) 
err_R_n = (A_hat-A_hat_n_R);
delta_R_n_i2 = mean(err_R_n.^2,2);

%%% Rank revealing QR bi-fidelity - with rank N_hi

% Normalize matrices for bound. 
A_samp = (A(:,1:n_est))./norm(A(:,1:n_est), 'fro');
B_samp = (B)./norm(B, 'fro');

% Obtain column skeleton of P
[P_s,ix] = matrixIDvR(B_samp,N_hi);

% Sample to compute bound - use ix + additional to reach n
rand_sample = [ix, getfield(setxor(ix,1:n), {1:n-numel(ix)})];

%%% rho_k - rank revealing qr

% Error bound inputs
normC = norm(P_s);
sb = svd(B_samp); 
err_Bbar = norm(B_samp-B_samp(:,ix)*P_s); 

% calculate true error to compare to A_hat and rho_k
err_Abar = norm(A_samp-A_samp(:,ix)*P_s)*norm(A(:,1:n_est), 'fro'); 
R = rank(A_samp-A_samp(:,ix)*P_s); 

% Optional: compare err_Abar to err_Ahat (MID to SMR)
% err_Ahat = norm(A_hat-A(:,1:n_est));
% [err_Ahat, err_Abar]

% Pre-process operation
B_R = B_samp(:,rand_sample);
A_R = A_samp(:,rand_sample);

% Scale MID estimate by norm(A(:,sample), 'fro')

% Compute epsilon tau
[~, rho_k,~, ~,~] = ...
    mat_id_error_est_one_normal(B_R, A_R, normC, err_Bbar, sb, N_hi, n);
% not computing bound correctly - rho_k_1 should be less than err_Abar. 

% Calculate eqn (16) / Theorem 1.
rho_k = rho_k*norm(A(:,1:n_est), 'fro'); 

MID_efficacy = rho_k/err_Abar;

%%% zeta_i
%%% Option 1 estimate zeta_i from (22)
zeta_i_1 = err_bi_mean/term_1*n_est/rho_k^2;  

%%% Option 2 estimate zeta_i from (24)
zeta_i_2 = delta_R_N_i2*n_est/rho_k^2;  

% n_est is no longer right? 
zeta_iN_2 = delta_R_n_i2*n_est/rho_k^2;  

%%% \bar{zeta}_N
%%% Option 1 estimate zeta_i from (23)
zeta_N_1 = err_bi_sum/term_1*n_est/R/rho_k^2; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound via moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eqn (33) from BiFiPCE_JH
V = (A_hat(:,sample)-A(:,sample)).^2; 
W = sum(V); 

% calculate moments
mu_V = mean(V,2); %1e-11
sigma_V = mean(abs(V-mu_V),2); % possibly meaningless - machine precision (-50)
rho_V = mean(abs(V-mu_V),2).^3; % -e-76 and some negatives... 

mu_W = mean(W); 
sigma_W =  mean(abs(W-mu_W),2);
rho_W =  mean(abs(W-mu_W),2).^3;

phi_t = normcdf(t); 

p_33 = phi_t - C.*rho_V./(sigma_V.^3*sqrt(N_hi));
bound_34 = term_1.*(mu_V+t.*sigma_V/sqrt(N_hi)); 

p_35 = phi_t - C.*rho_W./(sigma_W.^3*sqrt(N_hi));
bound_36 = term_1*(mu_W + t*sigma_W/sqrt(N_hi)); 

eff_36 = sqrt(bound_36/err_bi_sum); 


