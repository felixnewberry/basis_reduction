% C_tilde as function of n. 

clear all
close all
clc

tic 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend

% Paper size
size_1 = [0,0,575,445]; 
size_2 = [0,0,1150,445]; 

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


% airfoil? 
% LDC? 
% what are respective n, r and N values? 




%%% Airfoil 
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

n_points = 200; 
x_int = linspace(-1,1,n_points); 

n_est = 100; 
% n_est = 500; 

xi_low = xi_ref(1:n_est,:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Params
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RLS update
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Want to model eqn (5) where n changes. %
% Don't care about coefficients - that's the solution. 
% Right now We're interested in eta_n which corresponds to A in RLS
% C = (A^T A)^-1 
% C_tilde^-1 = C^01 + WWT 
% All I really need is eta, as additional samples are appended. That
% shouldn't be too hard shoudl it? 
% A little bit hard - eta is a function of alpha psi - I have to do some
% math to pin it down. 

% Just follow through bi-fid algorithim, but change n and compute new eta 
% For each eta, compunte C and plot. 
% Plot the values of C or norm of C? Decide once I have C. 

% can i set up eta to more clearly append new measurements? Inspect n and
% n+1

r = 5; 
R = r+10;
n = 5; 

n_step = 20;
% N_hi_vec = [5,6,7,8,9,10:n_step:n_est]; 
N_hi_vec = [n,10,20:n_step:500]; 
% N_hi_vec = 5; 

norm_C_vec = zeros(length(N_hi_vec),1); 
norm_C_diff_vec = zeros(length(N_hi_vec),1); 

norm_CW_vec = zeros(length(N_hi_vec)-1,1); 
% norm_CW_diff_vec = zeros(length(N_hi_vec)-1,1); 

for i_N_hi = 1:length(N_hi_vec)
    % each addition of a sample appends 1 new value to psi_bi_n
    [psi_bi_n, psi_bi] = br_eta_n(A_inf, N_hi_vec(i_N_hi), psi_ref, c_low, sigma, r, n_est);
    
    W_n = psi_bi_n(n+1:end,2:end)'; 
    
    if i_N_hi == 1
        A_N = psi_bi(:,2:end); 
        C_N = inv(A_N'*A_N);
        
        W_N = psi_bi(n+1:end,2:end)'; 
%         C_N_W_n = (A_N'*A_N)\W_N;
    end
    A_n = psi_bi_n(:,2:end); 
    C_n = inv(A_n'*A_n);
    
    norm_C_vec(i_N_hi) = norm(C_n); 
    norm_C_diff_vec(i_N_hi) = norm(C_n - C_N); 
    
    if i_N_hi >= 2
        norm_CW_vec(i_N_hi-1) = norm((A_n'*A_n)\W_n); 
%         norm_CW_diff_vec(i_N_hi-1) = norm((A_n'*A_n)\W_n - C_N_W_n); 
    end
end

% I could also store all the values... 

% make plot pretty - add title

% Is this converging? Make plot pretty, and test 200. Then LDC as well.
% Should I look at values too? for r = 5, 15 values. how about LDC? 
figure
semilogy(N_hi_vec, norm_C_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert \tilde{C} \Vert$ ', 'interpreter', 'latex', 'fontsize', FS)
axis tight
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
title('Airfoil','interpreter', 'latex', 'fontsize', FS_leg)

figure
p1 = semilogy(N_hi_vec, norm_C_diff_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert C_n - \tilde{C} \Vert$ ', 'interpreter', 'latex', 'fontsize', FS)
axis tight
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
title('Airfoil','interpreter', 'latex', 'fontsize', FS_leg)

figure
semilogy(N_hi_vec(2:end), norm_CW_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert \tilde{C} W\Vert$ ', 'interpreter', 'latex', 'fontsize', FS)
axis tight
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
title('Airfoil','interpreter', 'latex', 'fontsize', FS_leg)

% figure
% p1 = semilogy(N_hi_vec(2:end), norm_CW_diff_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
% axis tight
% xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$\Vert C_n - \tilde{C} \Vert$ ', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)

