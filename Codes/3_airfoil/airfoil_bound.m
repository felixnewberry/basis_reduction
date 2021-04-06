% Airfoil 

% Bound break down

clear all
close all
clc

tic 

save_on = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot Settings                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend

size_1 = [0,0,445,345]; 
size_2 = [0,0,1340,515]; 

size_square = [0,0,445,445]; 
size_large = [0,0,668,518]; 

FS = 28;    % Font size axis
FS_axis = 18; 
LW_axis = 1; 

% Colors
c1 = [0, 0.4470, 0.7410];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [0.9290, 0.6940, 0.1250]; 
c4 = [0.4940, 0.1840, 0.5560];
c5 = [0.4660, 0.6740, 0.1880]; 
c6 = [0.3010, 0.7450, 0.9330]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Airfoil 
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
x_l = linspace(-1,1,n_points); 

n_est = 500; 
% n_est = 200; 

xi_low = xi_low(1:n_est,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 5;                          % PCE order
d = 6;                          % Stochastic dimension

% tolerance on residual used in spgl1 solver
sigma = .001;

pc_solver = 2;  %0 is LS, 1 is mmv and 2 is individual spg
pc_val = 0; 

n_sim = length(xi_ref); 

t_setup = toc; 
fprintf('Setup with time : %d s.\n', t_setup);
timerval = tic;

% Statistical error bound params - review C
C = 0.4748; 
% Choose t and corresponding phi_t (cumulative distribution for standard
% normal. 
t = 2; 

n_reps = 30; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uc = u_low'; Uf = u_ref'; 
B = Uc(:,1:n_est)/norm(Uc(:,1:n_est),'fro');
A = Uf(:,1:n_est)/norm(Uf(:,1:n_est),'fro');

A_inf = Uf/norm(Uf(:,1:n_est),'fro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% reference solution - psi_ref is used... 
% [c_ref, psi_ref] = my_pce(xi_ref, p, A', sigma, pc_solver, pc_val); 

[c_ref, psi_ref] = my_pce(xi_ref, p, A_inf', sigma, pc_solver, pc_val); 
t_ref = toc(timerval); 
fprintf('Reference solution : %d s.\n', t_ref);
timerval = tic; 

%%% Low-fidelity solution
[c_low, psi_low] = my_pce(xi_low, p, B', sigma, pc_solver, pc_val); 
t_low = toc(timerval); 
% t1 = toc
% t_ref
fprintf('Low fidelity solution : %d s.\n', t_low);


%%% Low_fidelity 
% u_low_int = interp1(x_l, B, x_h);  
err_low = norm(B-A)/norm(A); 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Efficacy loop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Where is the parfor stuff? Maybe I didn't do that for these... which is
% % silly. Try run the LDC now and see how it goes? 
% N_hi_vec = 3:50; 
% r_vec = 3:20; 
% 
% % N_hi_vec = 30; 
% % r_vec = 8; 
% 
% %r_vec = [3, 20]; 
% %N_hi_vec = [3,50];
% %n_reps = 2;
% 
% length_n = length(N_hi_vec);
% length_r = length(r_vec); 
% 
% n_r_results{n_reps} = [];
% 
% % Repeat
% parfor i_rep = 1:n_reps
% % for i_rep = 1:n_reps
%     i_rep
%     n_r_results{i_rep}.efficacy = zeros(length_r,length_n); 
%     n_r_results{i_rep}.prob = zeros(length_r,length_n); 
% 
%     for i_n = 1:length_n
%         for i_r = 1:length_r
% 
%         [n_r_results{i_rep}.efficacy(i_r, i_n), n_r_results{i_rep}.prob(i_r, i_n)] = ...
%     br_bound_practical(A_inf, N_hi_vec(i_n), psi_ref, c_low, sigma, r_vec(i_r), n_est, t, C);
% 
%         end
%     end
% end
% 1;
% 
% % Calculate mean over n_reps
% stat_struct = cat(3,n_r_results{:});
% efficacy_mat = mean(cat(3,stat_struct.efficacy),3);
% prob_mat = mean(cat(3,stat_struct.prob),3);
% 
% % efficacy_mat = mean_ep_tau_bound./mean_bi_err; 
% 
% if save_on == 1
%     save('Results/Airfoil_efficacy_1', 'r_vec', 'efficacy_mat', 'prob_mat', 'N_hi_vec')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single point:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = 8; 
N_hi = 30; 
R = N_hi+10;

% err_bi_mean_rep     = zeros(n_points,n_reps);
% err_bi_sum_rep      = zeros(1, n_reps); 
% mu_rep              = zeros(1, n_reps);
% rho_k_rep           = zeros(1, n_reps);
% zeta_vec_rep        = zeros(n_points,n_reps);
% zeta_num_rep        = zeros(n_points,n_reps);
% zeta_den_rep        = zeros(n_points,n_reps);
% theta_vec_rep       = zeros(n_points,n_reps);
% MID_efficacy_rep    = zeros(1, n_reps);
% eff_vec_rep         = zeros(5, n_reps);
% factor_vec_rep      = zeros(4, n_reps);
% bound_34_rep        = zeros(n_points,n_reps);
% bound_20_rep        = zeros(n_points,n_reps);
% p_33_rep            = zeros(n_points,n_reps);
% p_35_rep            = zeros(1,n_reps);
% R_rep               = zeros(1,n_reps);

err_bi_mean_rep     = zeros(n_points,n_reps);
err_bi_sum_rep      = zeros(1, n_reps); 
mu_rep              = zeros(1, n_reps);
rho_k_rep           = zeros(1, n_reps);
zeta_i_1_rep        = zeros(n_points,n_reps);
zeta_i_2_rep        = zeros(n_points,n_reps);
zeta_N_1_rep        = zeros(1,n_reps);
% theta_vec_rep       = zeros(n_points,n_reps);
MID_efficacy_rep    = zeros(1, n_reps);
eff_36_rep         = zeros(1, n_reps);
% factor_vec_rep      = zeros(4, n_reps);
bound_34_rep        = zeros(n_points,n_reps);
% bound_20_rep        = zeros(n_points,n_reps);
p_33_rep            = zeros(n_points,n_reps);
p_35_rep            = zeros(1,n_reps);
R_rep               = zeros(1,n_reps);

for i_rep = 1:n_reps
    
    i_rep
    [err_bi_mean_rep(:,i_rep), err_bi_sum_rep(:,i_rep), mu_rep(:,i_rep),...
        rho_k_rep(:,i_rep), zeta_i_1_rep(:,i_rep), zeta_i_2_rep(:,i_rep), ...
        zeta_N_1_rep(:,i_rep), MID_efficacy_rep(:,i_rep),...
        eff_36_rep(:,i_rep), bound_34_rep(:,i_rep), ...
        p_33_rep(:,i_rep), p_35_rep(:,i_rep), R_rep(:,i_rep)]...
    =  br_bound_complete(B, A_inf, N_hi, R, psi_ref, c_low, sigma, r, n_est, t, C);
end

1; 

err_bi_mean     = mean(err_bi_mean_rep,2); 
err_bi_sum      = mean(err_bi_sum_rep);
mu              = mean(mu_rep);
rho_k           = mean(rho_k_rep);
zeta_i_1        = mean(zeta_i_1_rep,2); 
zeta_i_2        = mean(zeta_i_2_rep,2); 
zeta_N        = mean(zeta_N_1_rep,2); 
% theta_vec       = mean(theta_vec_rep,2); 
MID_efficacy    = mean(MID_efficacy_rep);
eff_36         = mean(eff_36_rep,2); 
% factor_vec      = mean(factor_vec_rep,2); 
bound_34        = mean(bound_34_rep,2); 
% bound_20        = mean(bound_20_rep,2); 
p_33            = mean(p_33_rep,2); 
p_35            = mean(p_35_rep);
R_rank          = mean(R_rep);

eff_36

if save_on == 1
    load('Results/Airfoil_efficacy_1');
    save('Results/Airfoil_bound_results_theta_est', 'err_bi_mean',...
    'N_hi_vec', 'r', 'R', 'N_hi', 'n_reps', 'x_l', ...
    'err_bi_sum', 'mu', 'rho_k', 'zeta_i_1', 'zeta_i_2', 'zeta_N',...
    'MID_efficacy', 'eff_36', 'bound_34', ...
    'p_33', 'p_35', 'R_rank');
end
   
