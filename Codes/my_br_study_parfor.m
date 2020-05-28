function [bi_stats, mean_lam_hi, lam_ref, lam_low]...
    = my_br_study(r,N_hi, n_reps, u_ref, xi_ref, psi_ref, sigma, c_low, c_ref)

%Loop over r and N_hi and construct bi-fidelty estimate and reference
%solution

%%% Inputs: 
% r         % KL truncation terms (vector to step through)
% N_hi      % Number of high fidelity samples (vector to step through)
% n_reps    % Number of repetitions
% u_ref     % Reference solution
% xi_ref    % Reference stochastic inputs

%%% Outputs: 

%%% define lengths for use in parfor
length_r = length(r);
length_N_hi = length(N_hi); 

% n_points = size(u_ref,2); 
n_samps = size(u_ref,1); 

1; 

% reference statistics
mean_ref = c_ref(:,1); 
var_ref = sum(c_ref(:,2:end).^2,2); 

%%% Stats that are independent of r and N_hi 

% Nodal covariance of the reference model
cov_ref = c_ref(:, 2:end)*c_ref(:,2:end)';
% Eigenvalue decomposition of nodal covariance
[~,lam_ref] = eigs(cov_ref);
lam_ref = diag(lam_ref);

% Nodal covariance of the low fidelity model
cov_low = c_low(:, 2:end)*c_low(:,2:end)';
% Eigenvalue decomposition of nodal covariance
[~,lam_low] = eigs(cov_low);
lam_low = diag(lam_low);

% %%% Initialize loop variables: 
% mean_mean_bi_err = zeros(length(r),length(N_hi));
% mean_mean_hi_err = zeros(length(r),length(N_hi)); 
% 
% var_mean_bi_err = zeros(length(r), length(N_hi)); 
% var_mean_hi_err = zeros(length(r), length(N_hi)); 
% 
% mean_var_bi_err = zeros(length(r), length(N_hi)); 
% mean_var_hi_err = zeros(length(r), length(N_hi)); 
% 
% var_var_bi_err = zeros(length(r), length(N_hi)); 
% var_var_hi_err = zeros(length(r), length(N_hi)); 
% 
% mean_lam_hi = cell(length(r), length(N_hi)); 
% mean_lam_ref = cell(length(r), length(N_hi)); 
% mean_lam_low = cell(length(r), length(N_hi)); 

disp('about to start r and N loop') 

1; 

% What is the data I want to arrive at? 

%%% Data
N_r_stats{n_reps} = [];
        
parfor i_rep = 1:n_reps    %index number of repetitions
    
%     % statistics        
%     N_r_stats.mean_bi_n = zeros(n_points, n_reps); 
%     N_r_stats.var_bi_n = zeros(n_points, n_reps); 
%     N_r_stats.mean_hi_n = zeros(n_points, n_reps); 
%     N_r_stats.var_hi_n = zeros(n_points, n_reps);         
% 
%     % eigenvalues
%     N_r_stats.lam_ref = zeros(r(i_rank), n_reps); 
%     N_r_stats.lam_hi = zeros(r(i_rank), n_reps); 
%     N_r_stats.lam_low = zeros(r(i_rank), n_reps); 
% 
%     % errors, hi and bi-fidelity
%     N_r_stats.mean_bi_err = zeros(1, n_reps); 
%     N_r_stats.mean_hi_err = zeros(1, n_reps); 
%     N_r_stats.var_bi_err = zeros(1, n_reps); 
%     N_r_stats.var_hi_err = zeros(1, n_reps); 

    % errors, hi and bi-fidelity
    N_r_stats{i_rep}.mean_bi_err = zeros(length_r, length_N_hi); 
    N_r_stats{i_rep}.mean_hi_err = zeros(length_r, length_N_hi); 
    N_r_stats{i_rep}.var_bi_err = zeros(length_r, length_N_hi); 
    N_r_stats{i_rep}.var_hi_err = zeros(length_r, length_N_hi); 
     
	% eigenvalues
    N_r_stats{i_rep}.lam_hi = cell(length_r, length_N_hi); 
    
%     
    for i_rank = 1:length_r        % index range of rank of bi-fid
        for i_hi = 1:length_N_hi   % index range of number of high fidelity simulations

            sample = datasample(1:n_samps, N_hi(i_hi), 'Replace', false); 
%             sample = 1:N_hi; 

            %%% High fidelity model - limited samples N_hi
            xi_hi = xi_ref(sample,:); 
            u_hi = u_ref(sample,:); 
            psi_hi= psi_ref(sample,:); 
            
            opts = spgSetParms('iterations',10000,'verbosity',0);
            
%             Hi fid PCE solved for via ell_1 minimization
            c_hi = spg_mmv(psi_hi,u_hi,sigma*norm(u_hi),opts);            
%             c_hi = psi_hi\u_hi; 
            c_hi = c_hi';
            
            mean_hi_n = c_hi(:,1);
            var_hi_n = (sum(c_hi(:,2:end).^2,2));
%             std_hi = sqrt((sum(c_hi(:,2:end).^2,2))); 
            
            %%% Bi-fidelity - many L, some H
            [mean_bi_n,var_bi_n,c_bi1,...
                lam_hi, psi_bi_low,alpha2]= ... 
                BR_FN(u_hi,psi_hi,c_hi,c_low,r(i_rank),sigma); 
            
%             %%% Bi-fidelity - many L, many H
%             [~,~,~,~,~,~, psi_bi_hi,~]= ...
%                 BR_FN(u_hi,psi_hi,c_hi,c_hi,c_ref,r(i_rank),sigma);
            
            N_r_stats{i_rep}.mean_bi_err(i_rank, i_hi) = norm(mean_bi_n - mean_ref)/norm(mean_ref); 
            N_r_stats{i_rep}.mean_hi_err(i_rank, i_hi) = norm(mean_hi_n - mean_ref)/norm(mean_ref); 
            N_r_stats{i_rep}.var_bi_err(i_rank, i_hi) = norm(var_bi_n- var_ref)/norm(var_ref); 
            N_r_stats{i_rep}.var_hi_err(i_rank, i_hi) = norm(var_hi_n - var_ref)/norm(var_ref);      
            
            
            
%             std_bi_err = norm(sqrt(var_bi_n(:,i_rep))-sqrt(var_ref))/norm(sqrt(var_ref));
            
%             fprintf('repetition : %d of %d \n', i_rep, n_reps);
%             fprintf('Timer : %d s.\n', toc);
            
            % eigenvalues for high
%             N_r_stats{i_rep}.lam_hi{i_rank, i_hi} = num2cell(lam_hi);
            N_r_stats{i_rep}.lam_hi{i_rank, i_hi} = lam_hi;
            
       end
        

    end
end

1; 




% Calculate mean over n_reps 
stat_struct = cat(3,N_r_stats{:});

mean_mean_bi_err = mean(cat(3,stat_struct.mean_bi_err),3);
mean_mean_hi_err = mean(cat(3,stat_struct.mean_hi_err),3);

var_mean_bi_err = var(cat(3,stat_struct.mean_bi_err),0,3);
var_mean_hi_err = var(cat(3,stat_struct.mean_hi_err),0,3);

mean_var_bi_err = mean(cat(3,stat_struct.var_bi_err),3);
mean_var_hi_err = mean(cat(3,stat_struct.var_hi_err),3);

var_var_bi_err = var(cat(3,stat_struct.var_bi_err),0,3);
var_var_hi_err = var(cat(3,stat_struct.var_hi_err),0,3);

% eigenvalues 
mean_lam_hi = cell(length_r, length_N_hi);

lam_struct = cat(3,stat_struct.lam_hi);
for i_r = 1:length_r
    for i_hi = 1:length_N_hi
        mean_lam_hi{i_r, i_hi} = mean(cat(2,lam_struct{i_r,i_hi,:}),2);
    end
end

% Group statistic results. Lambdas are seperate
bi_stats = cat(3,mean_mean_bi_err , mean_mean_hi_err, var_mean_bi_err, ...
    var_mean_hi_err, mean_var_bi_err , mean_var_hi_err ,...
    var_var_bi_err, var_var_hi_err);
end

