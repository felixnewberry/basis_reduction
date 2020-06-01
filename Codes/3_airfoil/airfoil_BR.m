clear all
close all
clc

%%% Airfoil 

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Choose grid points to look at, 128 total.
% gridpt = [1:128];
% 
% load('x_locations')
% x_l = x_locations; 
% x_h = x_locations;
% 
% load lowFiResults_6d2.mat;
% u_low = lowFiResults(gridpt,:)'; % Number of samples x spatial/temporal dimension of problem
% 
% load highFiResults_6d2.mat;
% u_ref = highFiResults(gridpt,:)';
% 
% % load random variables used to generate u
% load Inputs_6d2.mat;
% 
% %     normalize to be U~[-1,1] and re-label as xi.
% % xi_ref should be of dimension: Number of reference samples x stochastic
% % dimension d
% xi_ref =[(1/3)*alpha',(4*c-4)',(50*m-1)',(10*(mach-.2)-2)',(10*p -4)',((t - .1)*50 - 1)'];
% xi_low = xi_ref; % number of low fidelity samples X stochastic dimension d
% 
% 
% %%% Interpolate data to be evenly distributed
% % Adjust coordinates to go from -1 to 1. 
% % clockwise from trailing edge: positive is pressure surface (lower)
% 
% x_order = [1-x_h(65:end);x_h(1:64)-1];
% u_low_order = [u_low(:,65:end), u_low(:,1:64)];
% u_ref_order = [u_ref(:,65:end), u_ref(:,1:64)];
% 
% % u_low and u_ref
% n_points = 100; 
% n_sim = length(u_ref); 
% 
% x_int = linspace(-1,1,n_points); 
% u_low = zeros(n_sim, n_points); 
% u_ref = zeros(n_sim, n_points); 
% for i_sim = 1:n_sim
%     u_low(i_sim,:) = interp1(x_order, u_low_order(i_sim,:), x_int);
%     u_ref(i_sim,:) = interp1(x_order, u_ref_order(i_sim,:), x_int);
% end

load('airfoil_data_2/Uf.mat', 'Uf')
u_ref = Uf'; 

load('airfoil_data_2/cp_results_lims_pm5.mat')
u_low = cp_results(:,:,1)';

params = importdata(strcat('airfoil_data_2', '/setting_param_file.dat'), ',' , 1); 
params = params.data; 

% Standardize params to be -1 to 1. 

% m p, t and alpha - four parameters
m_lim = [0.032, 0.048]; % 0.4
p_lim = [0.32, 0.48]; 
t_nom = [0.096, 0.144]; 
alpha_nom = [0, 6];

lim_mat = [m_lim; p_lim; t_nom; alpha_nom]';


xi_ref = 2*(params(:,1:4)-lim_mat(1,:))./(lim_mat(2,:)-lim_mat(1,:))-1; 
xi_low = xi_ref; 

n_points = 200; 
x_int = linspace(-1,1,n_points); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 5;                          % PCE order
d = 6;                          % Stochastic dimension

N_hi = [10 30 50 70 90 110];      % Number high-fidelity samples
% N_hi = [30]; 

r = [3 8 10];               % KL order
% r = 3; 

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
var_low_err = norm(var_low - var_ref)/norm(var_ref); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity and High fidelity r and N_hi study 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[bi_stats, mean_lam_hi, mean_lam_ref, mean_lam_low]...
    = my_br_study(r, N_hi, n_reps, u_ref, xi_ref, psi_ref, sigma, c_low, c_ref, pc_solver);
% 
% % % Save results: 
% save('Results/Airfoil_results_int','bi_stats', 'mean_lam_hi', 'mean_lam_ref', ...
%     'mean_lam_low','N_hi',...
%     'var_low_err','mean_low_err', 'r')
% 
% save('Results/Airfoil_results_spg','bi_stats', 'mean_lam_hi', 'mean_lam_ref', ...
%     'mean_lam_low','N_hi',...
%     'var_low_err','mean_low_err', 'r')
% 
% save('Results/Airfoil_results_spg_int_2','bi_stats', 'mean_lam_hi', 'mean_lam_ref', ...
%     'mean_lam_low','N_hi',...
%     'var_low_err','mean_low_err', 'r')
