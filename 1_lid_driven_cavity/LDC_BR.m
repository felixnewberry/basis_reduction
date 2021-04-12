clear all
close all
clc

%%% Lid driven cavity

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('u_64_f_2.mat')
highFiResults = u_matrix_0'; 

Uc_nom_u = load('Nom_u_mid.mat', 'Uc','Ub','sb');
lowFiResults = Uc_nom_u.Uc; 

load 'x_64.mat'
x_h = x_64(:,1); 
x_l = x_h; 

load('xi_mat_2.mat')
xi_ref = xi_2; 
xi_low = xi_ref; 

% Grid points of interest
gridpt_l = 1:size(lowFiResults(:,1)); 
gridpt_h = 1:size(highFiResults(:,1)); 

% u should be number of samples x spatial/temporal dimension of the
% problem. size: n_samples x n_gridpoints
u_ref = highFiResults(gridpt_h,:)';
u_low = lowFiResults(gridpt_l,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 4;                          % PCE order
d = 2;                          % Stochastic dimension

N_hi = [10, 15, 20, 25];      % Number high-fidelity samples
% N_hi = [5]; 

r = [2 4 7];                  % KL order
% r = 4; 

% tolerance on residual used in spgl1 solver
sigma = .001;

% Number of repetitions
n_reps = 100; 

pc_solver = 2;  %0 is LS, 1 is mmv and 2 is individual spg

t_setup = toc; 
fprintf('Setup with time : %d s.\n', t_setup);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don't compute validation error. (Used for PC calibration)
pc_val = 0; 

%%% reference solution
[c_ref, psi_ref, ~] = my_pce(xi_ref, p, u_ref, sigma, pc_solver, pc_val); 

t_ref = toc - t_setup; 
fprintf('Reference solution : %d s.\n', t_ref);


%%% Low-fidelity solution

[c_low, psi_low, ~] = my_pce(xi_low, p, u_low, sigma, pc_solver, pc_val); 

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

% mean_low_int = interp1q(x_l', mean_low, x_h');
mean_low_err = norm(mean_low - mean_ref)/norm(mean_ref); 

% var_low_int = interp1q(x_l', var_low, x_h');
var_low_err = norm(mean_low - var_ref)/norm(var_ref); 

1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity r and N_hi study 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
1;

[bi_stats, mean_lam_hi, mean_lam_ref, mean_lam_low]...
    = my_br_study(r, N_hi, n_reps, u_ref, xi_ref, psi_ref, sigma, c_low, c_ref, pc_solver);


save('Results/LDC_results','bi_stats', 'mean_lam_hi', 'mean_lam_ref', ...
    'mean_lam_low','N_hi',...
    'var_low_err','mean_low_err', 'r')
