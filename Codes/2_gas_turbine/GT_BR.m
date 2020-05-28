clear all
close all
clc

%%% Gas Turbine
% Have to verify that low_20 is right not low_50... 

t_start = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Chose QoI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

QoI = 0; % u mid
% QoI = 1; % cylinder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% two options for now, 1/21 cost (assembledRunMid_40_2) or 1/51 (assembledRunMid_40) cost. 
% Checking using spg on each point for cylinder rather than mmv. Doing this
% for the reference and low, but not the hi within the bi-fidelity. Maybe
% it should  be there too. 
% Test each. 
% First with cylinder. 
    
if QoI == 0
    lowFiResults = importdata('assembledRunMid_40_1');   
    % I think _2 may have been used in earlier results?? 
%     lowFiResults = importdata('assembledRunMid_40_2'); 
    lowFiResults = lowFiResults';

    highFiResults = importdata('assembledRunMid_110_high'); 
    highFiResults = highFiResults';

    x_h = importdata('assembledCoordsMid_110_high');
    [~,idx_h] = sort(x_h(2,:));
    x_h = x_h(2,idx_h); 
    highFiResults = highFiResults(idx_h,:); 
  
    x_l = importdata('assembledCoordsMid_40_1');
%     x_l = importdata('assembledCoordsMid_40_2');
    [~,idx_l] = sort(x_l(2,:));
    x_l = x_l(2,idx_l);
    lowFiResults = lowFiResults(idx_l,:); 
    
elseif QoI == 1
    lowFiResults = importdata('assembledRunCylinder_40_1'); 
%     lowFiResults = importdata('assembledRunCylinder_40_2'); 
    lowFiResults = lowFiResults';
    tic

    highFiResults = importdata('assembledRunCylinder_110_high'); 
    highFiResults = highFiResults';
    
    x_h = importdata('assembledCoordsCylinder_110_high');
    x_h = atan2(x_h(2,:),-x_h(1,:));
    [~,idx_h] = unique(x_h);
    x_h = x_h(idx_h); 
    highFiResults = highFiResults(idx_h,:); 
  
    x_l = importdata('assembledCoordsCylinder_40_1');
%     x_l = importdata('assembledCoordsMid_40_2');
    x_l = atan2(x_l(2,:),-x_l(1,:));    
    [~,idx_l] = unique(x_l);
    x_l = x_l(idx_l);
    lowFiResults = lowFiResults(idx_l,:);  
end

y_samp = load('uniform_40_low.mat');
xi_low = y_samp.uniform_data;

y_samp = load('uniform_110_high.mat');
xi_ref = y_samp.uniform_data;
   

gridpt_l = 1:length(x_l); 
gridpt_h = 1:length(x_h); 

% gridpt_h = 1:4; 

u_ref = highFiResults(gridpt_h,:)';
u_low = lowFiResults(gridpt_l,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 6;                          % PCE order

N_hi = [10 30 50 70 90 110];    % Number high-fidelity samples
% N_hi = [30]; 

r = [3 8 10];                  % KL order
% r = 3; 

% tolerance on residual used in spgl1 solver
sigma = .007;

% Number of repetitions
n_reps = 100; 

pc_solver = 2;  %0 is LS, 1 is mmv and 2 is individual spg

t_setup = toc(t_start); 
fprintf('Setup with time : %d s.\n', t_setup);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% reference solution
[c_ref, psi_ref] = my_pce(xi_ref, p, u_ref, sigma, 2); 

t_ref = toc(t_start) - t_setup; 
fprintf('Reference solution : %d s.\n', t_ref);
% norm(psi_ref*c_ref'-u_ref)/norm(u_ref);

1; 

%%% Low-fidelity solution - changing this to spg may bungle results..???
%%% Need to make it consistent - fix soon. 

%%% More important that it is honest - mean is easy to estimate for the
%%% temperature across the midplane. No biggie if bi- and hi are bad. 
% Variance more useful. Different for cylinder. 
[c_low, psi_low] = my_pce(xi_low, p, u_low, sigma, 1); 

t_low = toc(t_start) - t_ref; 
fprintf('Low fidelity solution : %d s.\n', t_low);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - Test mmv vs spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [c_ref_mmv, psi_ref] = my_pce(xi_ref, p, u_ref, sigma, 1); 
% [c_ref_spg, psi_ref] = my_pce(xi_ref, p, u_ref, sigma, 2); 
% 
% [c_low_mmv, psi_low] = my_pce(xi_low, p, u_low, sigma, 1); 
% [c_low_spg, psi_low] = my_pce(xi_low, p, u_low, sigma, 2); 
% 
% mean_ref_mmv = c_ref_mmv(:,1); 
% mean_ref_spg = c_ref_spg(:,1); 
% 
% mean_low_mmv = c_low_mmv(:,1);
% mean_low_mmv_int = interp1q(x_l', mean_low_mmv, x_h');
% mean_low_spg = c_low_spg(:,1);
% mean_low_spg_int = interp1q(x_l', mean_low_spg, x_h');
% 
% % For the cylinder these appear to be the same. 
% % Likewise for the mid. 
% 
% % what about on mac?? 
% 1; 
% figure
% plot(mean_ref_spg,'k')
% hold on
% plot(mean_ref_mmv,'r')
% plot(mean_low_mmv_int,'b')
% plot(mean_low_spg_int,'g')
% 
% 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relative errors reference and low
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mean of reference solution
mean_ref = c_ref(:,1); 

% mean of low fidelity solution
mean_low = c_low(:,1);


var_ref = sum(c_ref(:,2:end).^2,2); 
std_ref = sqrt(sum(c_ref(:,2:end).^2,2)); 

% variance of low fidelity solution
var_low = sum(c_low(:,2:end).^2,2); 
std_low = sqrt(sum(c_low(:,2:end).^2,2)); 


mean_low_int = interp1q(x_l', mean_low, x_h');
mean_low_err = norm(mean_low_int - mean_ref)/norm(mean_ref); 

var_low_int = interp1q(x_l', var_low, x_h');
var_low_err = norm(var_low_int - var_ref)/norm(var_ref); 

mean_low_err


% solve low with spg_mmv, high with spg_bpdn for each point
% low: 0.9922
% If both are point specific, then: 
% 0.0041

% I think something is awry.. .in my_pce check this out. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity r and N_hi study 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Write a function

[bi_stats, mean_lam_hi, mean_lam_ref, mean_lam_low]...
    = my_br_study(r, N_hi, n_reps, u_ref, xi_ref, psi_ref, sigma, c_low, c_ref);

% Save results: 
if QoI == 0
    results_name = 'GT_u_mid_results';
elseif QoI == 1
    results_name = 'GT_cylinder_results'; 
end

mean_low_err
save(strcat('Results/',results_name),'bi_stats', 'mean_lam_hi', 'mean_lam_ref', ...
    'mean_lam_low','N_hi',...
    'var_low_err','mean_low_err', 'r')



% Have to adjust for time snapshot. and likely many other things. 


% GT n = 30, r = 3,8, 10 

error_mean_bi = bi_stats(1,1,1);
error_mean_hi = bi_stats(1,1,2);

error_var_bi = bi_stats(1,1,5);
error__var_hi = bi_stats(1,1,6);

% res = [mean_low_err, error_mean_hi, error_mean_bi; var_low_err, error__var_hi, error_var_bi]
% r= 3 n = 30

% I anticipate results for mid as 
% 2.7e-3, 3.1e-3, 2.3e-3
% 0.8e-1, 1.15e-1 1.15e-1

% I anticipate results for cylinder as 
% 1e-2, 1.3e-3, 1.05e-3
% 1e-1, 1.1e-1, 0.8e-1


% u mid using
%     0.0016    0.0130    0.0057
%     0.0813    0.5634    0.3606
   

% % % I think cylinder 50 (whatever that means, p = 6) is a goer. 

% I should print out the values that were plotted, then try match them - at least the low-fidelity error! 
% A worry that it doesn't already... 

% Test: 
% % % saves everything
% save('BR_FN_GT_Mid_low2_p6','mean_mean_bi_err' , 'mean_mean_hi_err' ,...
%     'var_mean_bi_err', 'var_mean_hi_err', ...
%     'mean_var_bi_err' , 'mean_var_hi_err' ,...
%     'var_var_bi_err', 'var_var_hi_err', ...
%     'mean_lam_hi', 'mean_lam_ref', 'mean_lam_low','N_hi',...
%     'var_low_err','mean_low_err','r', 'pol')

