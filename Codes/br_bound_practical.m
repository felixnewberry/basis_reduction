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
% err_bi_sum    - bi-fidelity error squared sum
% err_bi_vec    - bi-fidelity error for each point squared
% mu            - coherence
% rho_k         - bi-fid error rank revealing QR
% Y_Nh_vec      - Y for each point in M_H
% theta_vec     - theta for each point in M_H
% U_bar         - evaluation of norm in (27) - squared

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PC Bi-fidelity and High fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% rho_k is too large 
sample = datasample(1:n_est, N_hi, 'Replace', false); 
% sample = 1:N_hi; 

% sample = [ix, getfield(setxor(ix,1:N_hi), {1:N_hi-numel(ix)})];
% ensure sample uses ix from QR - and overlaps as much as possible with R 

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
%psi_bi_ref = psi_ref(:,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi = [ones(size(A(:,1:n_est), 2), 1) psi_bi]; 
% remove r+1 column before solving for coefficients.
% psi_bi = psi_bi(:,1:end); 
A_hat_n = (psi_bi*c_bi')'; 

% % If I calculated the norm over all n_est samples: 
% err_bi_vec = vecnorm(A_hat_n-A(:,1:n_est), 2, 2).^2;
% err_bi_sum = sum(err_bi_vec);

% calculate error for each sample and each point
err_bi_mat = (A_hat_n-A(:,1:n_est)).^2;
% error summed across all points - once for each sample. 
err_bi_sum = sum(err_bi_mat);

% err_bi_mean = mean(err_bi_mat,2); 

% To check that bi-fid estimate of just sample is similar to
% A_hat_n(:,sample)

% psi_bi2 = psi_ref(sample,2:end)*alpha2';
% psi_bi2 = [ones(size(A(:,sample), 2), 1) psi_bi2]; 
% A_hat_n2 = (psi_bi2*c_bi')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Error bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % estimate inf samples with inf samples
% [~, ~, c_bi_inf,~, ~, alpha2_inf]= BR_FN(A',psi_ref,c_hi,c_low,r,sigma); 
% psi_bi_inf = psi_ref(:,2:end)*alpha2_inf';
% psi_bi_inf = [ones(size(A(:,:), 2), 1) psi_bi_inf]; 
% A_hat_inf = (psi_bi_inf*c_bi_inf')'; 
% 
% %%% Coherence 
%  
% % psi_bi for all samples? No  - just for N_hi
% % A1 = psi_bi(:,2:end); 
% % mu1 = max(sum(abs(A1).^2,2));
A1 = psi_ref(sample,2:end)*alpha2';
mu = max(sum(abs(A1).^2,2)); 
% A1 is N_h x r not n_est x r

%%% X_bar

% evaluate theta from eqn (18)
% ratio of estimate with N_H samples compared to N_H=inf (or as close as we
% can get)
% Consider theta as necessary correction.

%  U_H changes size, we compare N_H to inf
% theta_vec = (vecnorm(A_hat_inf-A,2,2).^2./size(A,2))./(vecnorm(A_hat_n(:,sample)-A(:,sample),2,2).^2./N_hi);

% A_hat_inf is 65 x 2200 
% Check that A_hat_n(:,sample) is within machine precision of an
% estimate using just those samples - I think that's good enough. 
% norm(A_hat_n(:,sample) - A_hat_n2)

% comprised of theta_bar x Y

% eqn (21) to calculate Y
% err_bi = norm(A_hat_n-A(:,1:n_est));
% err_Bhat = norm(B-B(:,ix)*P_s); 

% Y is found from one instance (I think?)
% - i so calculate Y for each point and take average

% - the size of U_H is M_H by N_H... 

%%% Rank revealing QR bi-fidelity - estimate just N_h sample


% % Obtain column skeleton of P
% A_samp = (A(:,sample))./norm(A(:,sample), 'fro');
% B_samp = (B(:,sample))./norm(B(:,sample), 'fro');
% 
% [P_s,ix] = matrixIDvR(B_samp,r);
% 
% % rand_sample = [ix, getfield(setxor(ix,1:N), {1:n-numel(ix)})];
% % rand_sample = [ix, getfield(setxor(ix,1:n), {1:n-numel(ix)})];
% n_qr_bound = min(n, N_hi); % possibly this is an issue? % Don't sample more than estimate...?
% rand_sample = [ix, getfield(setxor(ix,1:n_qr_bound), {1:n_qr_bound-numel(ix)})];

% Have rand_sample include ix + some additional samples - not N_H? 
% Probably should be N_H - an improvement to make. 
% maybe random sample 

%%% rho_k - rank revealing qr to estimate just N_h samples = sample

% % Error bound inputs
% normC = norm(P_s);
% sb = svd(B_samp); 
% err_Bbar = norm(B_samp-B_samp(:,ix)*P_s); 
% 
% % calculate true error to compare to A_hat and rho_k
% err_Abar = norm(A_samp-A_samp(:,ix)*P_s)*norm(A(:,sample)); 
% 
% % compare to pc bi-fid: 
% err_Ahat = norm(A_hat_n(:,sample)-A(:,sample));


% Abar is smaller - but not by much ratio is err_Ahat/err_Abar = 1.4

% Subset of vectors for bi-fidelity error estimate 
% Ensure column skeleton selection are used + additional samples to reach
% total of n - done prior to pc bi-fid so that same ix samples are used
% in bi-fid. This seems more consistent.

% % Pre-process operation - normalize by frobenius norm or something like
% % that? 
% B_R = B_samp(:,rand_sample);
% A_R = A_samp(:,rand_sample);

% Now I need to scale things correctly: multiply by norm(A(:,sample), 'fro')

% % Compute epsilon tau
% [~, rho_k_1,~, ~,~] = ...
%     mat_id_error_est_one_normal(B_R, A_R, normC, err_Bbar, sb, N_hi, n);
% % not computing bound correctly - rho_k_1 should be less than err_Abar. 
% 
% % scale:
% rho_k = rho_k_1*norm(A(:,sample), 'fro');

% % Calculate Y from (21)
% Y_Nh_vec = (vecnorm(A(:,sample)-A_hat_n(:,sample),2,2).^2)./...
%     (vecnorm((A_samp-A_samp(:,ix)*P_s)*norm(A(:,sample), 'fro'), 2, 2).^2);
% Y_num = vecnorm(A(:,sample)-A_hat_n(:,sample), 2, 2); 
% Y_den = vecnorm((A_samp-A_samp(:,ix)*P_s)*norm(A(:,sample), 'fro'), 2, 2); 
% %vecnorm(Y_num, 2, 2)
% err_mat = [err_Abar, err_Ahat, rho_k]; 

% U_bar_vec = vecnorm((A_samp-A_samp(:,ix)*P_s)*norm(A(:,sample), 'fro'), 2, 2);
% U_hat_vec = vecnorm(A_hat_n(:,sample)-A(:,sample),2,2); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% New statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

term_1 = (1+4*mu/N_hi);

% Eqn (33) from BiFiPCE_JH
V = (A_hat_n(:,sample)-A(:,sample)).^2; 
W = sum(V); 

% calculate moments
% mu_V = mean(V,2); %1e-11
% sigma_V = mean(abs(V-mu_V),2); % possibly meaningless - machine precision (-50)
% rho_V = mean(abs(V-mu_V),2).^3; % -e-76 and some negatives... 

mu_W = mean(W); 
sigma_W =  mean(abs(W-mu_W),2);
rho_W =  mean(abs(W-mu_W),2).^3;

% For a given probability (39) and (41) 

phi_t = normcdf(t); 

% % eqn (39)
% p_39 = phi_t - C.*rho_V./(sigma_V.^3*sqrt(N_hi));
% bound_40 = term_1.*(mu_V+t.*sigma_V/sqrt(N_hi)); 

p_41 = phi_t - C.*rho_W./(sigma_W.^3*sqrt(N_hi));
bound = term_1*(mu_W + t*sigma_W/sqrt(N_hi)); 

% efficacy_40 = sqrt(bound_40./mean(err_bi_mat,2));
% efficacy_42 = sqrt(bound_42./mean(err_bi_sum));
% rho_W./(sigma_W.^3*sqrt(N_hi)) % why is this always 0.2236? 

efficacy = sqrt(bound./mean(err_bi_sum));
end

