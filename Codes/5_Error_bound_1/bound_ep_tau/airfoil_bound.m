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
x_int = linspace(-1,1,n_points); 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uc = u_low'; Uf = u_ref'; 
B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% reference solution
[c_ref, psi_ref] = my_pce(xi_ref, p, A', sigma, pc_solver, pc_val); 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over N_hi, or n to compute epsilon tau, or r for reduced basis... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This should be put in a n_reps parfor loop

% Keep N_hi constant - change number of samples used in estimate (n_vec)
% Also change truncation r. 



% n_vec = N_hi-10:N_hi+10;
% n_vec = 3:4:20; 
N_hi_vec = 3:50; 



%%% Check that q > 0 
% works for 33 for LDC - Bound is very loose at this point. 
% Doesn't seem useful in general... 

% N_hi_vec = 3:4; 
% r = 3; 
r_vec = 3:20; 
n_vec = r_vec+10; 

length_n = length(N_hi_vec);
length_r = length(r_vec); 

n_reps =  30;


n_r_results{n_reps} = [];

% Repeat
for i_rep = 1:n_reps
    
    n_r_results{i_rep}.err_bi_vec = zeros(length_r,length_n); 
    n_r_results{i_rep}.err_ep_tau_bound_vec = zeros(length_r,length_n); 
    
    for i_n = 1:length_n
        for i_r = 1:length_r
        [n_r_results{i_rep}.err_bi_vec(i_r, i_n), n_r_results{i_rep}.err_ep_tau_bound_vec(i_r, i_n), ~] = ...
            br_ep_tau_error(B, A, N_hi_vec(i_n), n_vec(i_r), psi_ref, psi_low, c_low, sigma, ...
            r_vec(i_r), p, xi_low, pc_solver); 
        end
    end

end

% Calculate mean over n_reps
stat_struct = cat(3,n_r_results{:});
mean_bi_err = mean(cat(3,stat_struct.err_bi_vec),3);
mean_ep_tau_bound = mean(cat(3,stat_struct.err_ep_tau_bound_vec),3);

efficacy = mean_ep_tau_bound./mean_bi_err; 

1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save('bound_ep_tau/airfoil_efficacy', 'n_vec', 'r_vec', 'efficacy', 'n_reps', 'N_hi_vec')

figure
h = pcolor(N_hi_vec, r_vec, efficacy);
set(h, 'EdgeColor', 'none');
axis tight
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('rank $r$', 'interpreter', 'latex', 'fontsize', FS)
colorbar
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Efficacy','interpreter', 'latex', 'fontsize', FS_leg)

% if save_on == 1
%     saveas(gcf,'bound_plots/Airfoil_bound_efficacy_1','png')
% end 


% appropriate range? Try figure this out... 
% For airfoil: 
% Have a run-save... 

% For LDC

% FOr GT - two QoI to do. Run overnight. Have airfoil and LDC for
% tomorrow's meeting. 
