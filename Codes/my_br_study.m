function [bi_stats, mean_lam_hi, mean_lam_ref, mean_lam_low]...
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

% reference statistics
mean_ref = c_ref(:,1); 
var_ref = sum(c_ref(:,2:end).^2,2); 

mean_mean_bi_err = zeros(length(r),length(N_hi));
mean_mean_hi_err = zeros(length(r),length(N_hi)); 

var_mean_bi_err = zeros(length(r), length(N_hi)); 
var_mean_hi_err = zeros(length(r), length(N_hi)); 

mean_var_bi_err = zeros(length(r), length(N_hi)); 
mean_var_hi_err = zeros(length(r), length(N_hi)); 

var_var_bi_err = zeros(length(r), length(N_hi)); 
var_var_hi_err = zeros(length(r), length(N_hi)); 

mean_lam_hi = cell(length(r), length(N_hi)); 
mean_lam_ref = cell(length(r), length(N_hi)); 
mean_lam_low = cell(length(r), length(N_hi)); 

disp('about to start r and N loop') 

for i_rank = 1:length(r)        % index range of rank of bi-fid
    for i_hi = 1:length(N_hi)   % index range of number of high fidelity simulations
   
        % statistics        
        mean_bi_n = zeros(size(u_ref,2), n_reps); 
        var_bi_n = zeros(size(u_ref, 2), n_reps); 
        mean_hi_n = zeros(size(u_ref, 2), n_reps); 
        var_hi_n = zeros(size(u_ref, 2), n_reps);         
        
        % eigenvalues
        lam_ref = zeros(r(i_rank), n_reps); 
        lam_hi = zeros(r(i_rank), n_reps); 
        lam_low = zeros(r(i_rank), n_reps); 
        
        % errors, hi and bi-fidelity
        mean_bi_err = zeros(1, n_reps); 
        mean_hi_err = zeros(1, n_reps); 
        var_bi_err = zeros(1, n_reps); 
        var_hi_err = zeros(1, n_reps); 
        
       for i_rep = 1:n_reps    %index number of repetitions
            
            sample = datasample(1:size(u_ref,1), N_hi(i_hi), 'Replace', false); 
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
            
            mean_hi_n(:,i_rep) = c_hi(:,1);
            var_hi_n(:,i_rep) = (sum(c_hi(:,2:end).^2,2));
            std_hi = sqrt((sum(c_hi(:,2:end).^2,2))); 
            
            %%% Bi-fidelity - many L, some H
            [mean_bi_n(:,i_rep),var_bi_n(:,i_rep),c_bi1,...
                lam_ref(:,i_rep),lam_hi(:,i_rep),lam_low(:,i_rep), psi_bi_low,alpha2]= ... 
                BR_FN(u_hi,psi_hi,c_hi,c_low,c_ref,r(i_rank),sigma); 
            
%             %%% Bi-fidelity - many L, many H
%             [~,~,~,~,~,~, psi_bi_hi,~]= ...
%                 BR_FN(u_hi,psi_hi,c_hi,c_hi,c_ref,r(i_rank),sigma);
            
            mean_bi_err(i_rep) = norm(mean_bi_n(:,i_rep) - mean_ref)/norm(mean_ref); 
            mean_hi_err(i_rep) = norm(mean_hi_n(:,i_rep) - mean_ref)/norm(mean_ref); 
            var_bi_err(i_rep) = norm(var_bi_n(:,i_rep)- var_ref)/norm(var_ref); 
            var_hi_err(i_rep) = norm(var_hi_n(:, i_rep) - var_ref)/norm(var_ref);
%             std_bi_err = norm(sqrt(var_bi_n(:,i_rep))-sqrt(var_ref))/norm(sqrt(var_ref));
            
%             fprintf('repetition : %d of %d \n', i_rep, n_reps);
%             fprintf('Timer : %d s.\n', toc);
%             1; 
            
       end
        
         % Take average of repetitions
        mean_mean_bi_err(i_rank,i_hi) = mean(mean_bi_err);
        mean_mean_hi_err(i_rank,i_hi) = mean(mean_hi_err);
        
        var_mean_bi_err(i_rank,i_hi) = var(mean_bi_err);
        var_mean_hi_err(i_rank,i_hi) = var(mean_hi_err);
        
        mean_var_bi_err(i_rank,i_hi) = mean(var_bi_err);
        mean_var_hi_err(i_rank,i_hi) = mean(var_hi_err);
        
        var_var_bi_err(i_rank,i_hi) = var(var_bi_err);
        var_var_hi_err(i_rank,i_hi) = var(var_hi_err);
     
%         eigenvalues
%         Note: saved as a cell because the number of eigenvalues changes 
%         w.r.t rank
        mean_lam_hi{i_rank,i_hi} = num2cell(mean(lam_hi,2)); 
        mean_lam_ref{i_rank,i_hi} = num2cell(mean(lam_ref,2)); 
        mean_lam_low{i_rank,i_hi} = num2cell(mean(lam_low,2));    
                
        fprintf('N_h  : %d \n', N_hi(i_hi));
        fprintf('r : %d s.\n', r(i_rank));
    end
end

% Group statistic results. Lambdas are seperate
bi_stats = cat(3,mean_mean_bi_err , mean_mean_hi_err, var_mean_bi_err, ...
    var_mean_hi_err, mean_var_bi_err , mean_var_hi_err ,...
    var_var_bi_err, var_var_hi_err);
end

