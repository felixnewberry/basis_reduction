clear all
close all
clc

%%% Lid driven cavity

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I should re-run this with 64x64 grid as high, and try 8x8 as well as 4x4
% as low. 

results_file_name = 'LDC_results'; 

load lowFiResults_LDC
load xi_low_LDC

load highFiResults_LDC.mat
load xi_high_LDC.mat

xi_ref = xi_high; 


% Grid points of interest
gridpt_l = 1:size(lowFiResults(:,1)); 
gridpt_h = 1:size(highFiResults(:,1)); 

% u should be number of samples x spatial/temporal dimension of the
% problem. size: n_samples x n_gridpoints
u_ref = highFiResults(gridpt_h,:)';
u_low = lowFiResults(gridpt_l,:)';

%%% Low and high coordinates
n_cell_l = 8;
x_l = linspace(0,1,n_cell_l+1);
x_l = 0.5*(cos(pi*(x_l-1)/2)+1);

n_cell_h = 32;
x_h = linspace(0,1,n_cell_h+1);
x_h = 0.5*(cos(pi*(x_h-1)/2)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p = 4;                          % PCE order
d = 2;                          % Stochastic dimension

N_hi = [10, 15, 20, 25];      % Number high-fidelity samples
N_hi = [10]; 

r = [1 3 6];                  % KL order
r = 3; 

% tolerance on residual used in spgl1 solver
sigma = .001;

% Number of repetitions
n_reps = 100; 

pc_solver = 1;  %0 is LS, 1 is mmv and 2 is individual spg

t_setup = toc; 
fprintf('Setup with time : %d s.\n', t_setup);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% reference solution
[c_ref, psi_ref] = my_pce(xi_ref, p, u_ref, sigma, pc_solver); 

t_ref = toc - t_setup; 
fprintf('Reference solution : %d s.\n', t_ref);


%%% Low-fidelity solution

[c_low, psi_low] = my_pce(xi_low, p, u_low, sigma, pc_solver); 

t_low = toc - t_ref; 
fprintf('Low fidelity solution : %d s.\n', t_low);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity r and N_hi study 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Write a function

[bi_stats, mean_lam_hi, mean_lam_ref, mean_lam_low]...
    = my_br_study(r, N_hi, n_reps, u_ref, xi_ref, psi_ref, sigma, c_low, c_ref);

% Save results: 
save('Results/LDC_results','bi_stats', 'mean_lam_hi', 'mean_lam_ref', ...
    'mean_lam_low','N_hi',...
    'var_low_err','mean_low_err', 'r')


% Have to adjust for time snapshot. and likely many other things. 
