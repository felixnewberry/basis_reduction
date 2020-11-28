function [psi_bi_n, psi_bi] = ...
    br_eta_n(A, N_hi, psi_ref, c_low, sigma, r, n_est)
% Compute bi-fidelity error and error bound quantities

% Inputs: 
% A - highfidelity samples - first n_est correspond to B - then additional
% N_hi - number of high-fidelity samples
% psi_ref - pc basis of complete high fidelity samples
% sigma - tolerance on residual in spgl1 solver
% c_low - low fidelity pc coefficients
% r - truncation of KL expansion 

% Outputs: 
% eta_bi - bi-fidelity basis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample = datasample(1:size(A,2), N_hi, 'Replace', false); 

sample = 1:N_hi; 

%%% High fidelity model - limited samples N_hi
% xi_hi = xi_ref(sample,:); 
u_hi = A(:,sample)'; 
psi_hi= psi_ref(sample,:); 

% opts = spgSetParms('iterations',10000,'verbosity',0);
% 
% %             Hi fid PCE solved for via ell_1 minimization
% c_hi = spg_mmv(psi_hi,u_hi,sigma*norm(u_hi),opts);            
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
                    
            
[~, ~, c_bi,~, psi_bi_n, alpha2]= BR_FN(u_hi,psi_hi,c_hi,c_low,r,sigma); 
            
% %%% Bi-fidelity error calculation -for estimate
% psi_bi = psi_ref(1:n_est,2:end)*alpha2';
% %psi_bi_ref = psi_ref(:,2:end)*alpha2';
% % New reduced basis including column of ones
% psi_bi = [ones(size(A(:,1:n_est), 2), 1) psi_bi]; 
% % remove r+1 column before solving for coefficients.
% psi_bi = psi_bi(:,1:end); 

%%% Bi-fidelity error calculation - for n tends to infty
psi_bi = psi_ref(:,2:end)*alpha2';
%psi_bi_ref = psi_ref(:,2:end)*alpha2';
% New reduced basis including column of ones
psi_bi = [ones(size(A(:,:), 2), 1) psi_bi]; 
% remove r+1 column before solving for coefficients.
psi_bi = psi_bi(:,1:end); 

%A_hat_n = (psi_bi*c_bi')'; 
%err_bi = norm(A_hat_n-A(:,1:n_est));

end

