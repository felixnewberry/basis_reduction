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
%%% Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = 5; 
R = r+10;

% sampling less then r doesn't make much sense... what to do?
N_hi_vec = [5,6,7,8,9,10:5:40]; 

n_rep = 10; 

err_bi = zeros(n_rep,length(N_hi_vec));
err_ep_tau_bound = zeros(1,length(N_hi_vec));
err_Bhat = zeros(1,length(N_hi_vec));
err_E = zeros(1,length(N_hi_vec));
err_T = zeros(1,length(N_hi_vec));
err_E_exact = zeros(1,length(N_hi_vec));
err_T_exact = zeros(1,length(N_hi_vec));
err_CS_pre = zeros(1,length(N_hi_vec));
err_CS_post = zeros(1,length(N_hi_vec));
err_bi_N = zeros(n_rep,length(N_hi_vec));



for i_rep =1:n_rep
    
for i_N_hi = 1:length(N_hi_vec)
    [err_bi(i_rep, i_N_hi), err_ep_tau_bound(i_N_hi), err_Bhat(i_N_hi),...
        err_E(i_N_hi), err_T(i_N_hi), Ir, err_E_exact(i_N_hi),...
        err_T_exact(i_N_hi), err_CS_pre(i_N_hi), err_CS_post(i_N_hi), err_bi_N(i_rep,i_N_hi)] = ...
                br_ep_tau_error_dev(B, A, N_hi_vec(i_N_hi), R, psi_ref, psi_low, c_low, sigma, ...
                r, p, xi_low, pc_solver);
end

end
% note that n samps in bound is fixed. Right now I'm changing n high
% sampes. ie R is fixed n is changing. 

figure
p1 = plot(N_hi_vec, err_bi,'-o','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(N_hi_vec, err_ep_tau_bound,'--x','Color',c2,...
    'LineWidth',LW,'MarkerSize',MS);
p3 = plot(N_hi_vec, 2*err_E_exact+err_T_exact.*err_Bhat,'-.s','Color',c3,...
    'LineWidth',LW,'MarkerSize',MS);
p4 = plot(N_hi_vec, err_E_exact+err_T_exact.*err_Bhat+err_CS_pre,':d','Color',c4,...
    'LineWidth',LW,'MarkerSize',MS);
p5 = plot(N_hi_vec, err_ep_tau_bound+err_bi_N,'-<','Color',c5,...
    'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1(1),p2, p3, p4,p5(1)],{'True', '$\epsilon ( \tau )$','E, T Exact', 'Pre CS, E, T Exact','Full Bound'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

% Interested in difference between H_hat_N and H_hat_n


figure
p5 = plot(N_hi_vec, err_bi_N,'-<','Color',c5,...
    'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
              