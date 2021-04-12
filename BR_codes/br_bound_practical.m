function [efficacy, p_41] = ...
    br_bound_practical(A, N_hi, psi_ref, c_low, sigma, r, n_est, t, C)
% Compute bi-fidelity error and error bound quantities

% Inputs: 
% B - low-fidelity samples - normalized by fro of first n_est. 
% A - highfidelity samples - first n_est correspond to B - then additional
% N_hi - number of high-fidelity samples
% n - number of samples used to compute bound
% psi_ref - pc basis of complete high fidelity samples
% sigma - tolerance on residual in spgl1 solver
% c_low - low fidelity pc coefficients
% r - truncation of KL expansion 

% Outputs: 
% efficacy 
% probability

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PC Bi-fidelity and High fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rho_k is too large 
sample = datasample(1:n_est, N_hi, 'Replace', false); 

%%% High fidelity model - limited samples N_hi
u_hi = A(:,sample)'; 
psi_hi= psi_ref(sample,:); 

weights = get_matrix_weights(psi_hi);
Psiweights = psi_hi*weights;    
norm_u_vec = vecnorm(u_hi,2); 

n_points = size(u_hi,2); 
P = size(psi_hi,2); 

c_data{n_points} = []; 
%parfor i_points=1:n_points
for i_points=1:n_points
    opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

    delta = sigma*norm_u_vec(i_points);
    c_data{i_points}.c_vec = weights*spg_bpdn(Psiweights,u_hi(:,i_points),delta,opts);
end

c_hi = zeros(n_points, P);
for i_points=1:n_points
    c_hi(i_points,:) = c_data{i_points}.c_vec; 
end


%%% High fidelity coeficients are used only to compute eigen values.    
[~, ~, c_bi,~, ~, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma); 
            
%%% Bi-fidelity error calculation
psi_bi = psi_ref(1:n_est,2:end)*alpha2';

% New reduced basis including column of ones
% remove r+1 column before solving for coefficients.
psi_bi = [ones(size(A(:,1:n_est), 2), 1) psi_bi(:,1:end-1)]; 

A_hat_n = (psi_bi*c_bi')'; 

% calculate error for each sample and each point
err_bi_mat = (A_hat_n-A(:,1:n_est)).^2;
% error summed across all points - once for each sample. 
err_bi_sum = sum(err_bi_mat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%% Coherence 
A1 = psi_ref(sample,2:end)*alpha2';
mu = max(sum(abs(A1).^2,2)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound via moments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term_1 = (1+4*mu/N_hi);

% Eqn (33) from BiFiPCE_JH
V = (A_hat_n(:,sample)-A(:,sample)).^2; 
W = sum(V); 

% calculate moments
mu_W = mean(W); 
sigma_W =  mean(abs(W-mu_W),2);
rho_W =  mean(abs(W-mu_W),2).^3;

phi_t = normcdf(t); 
 
p_41 = phi_t - C.*rho_W./(sigma_W.^3*sqrt(N_hi));
bound = term_1*(mu_W + t*sigma_W/sqrt(N_hi)); 

efficacy = sqrt(bound./mean(err_bi_sum));
end

