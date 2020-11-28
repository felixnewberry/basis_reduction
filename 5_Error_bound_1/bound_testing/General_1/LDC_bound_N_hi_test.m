% Test Error Bound

% Problem setup

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

% Lid driven cavity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res_name = 'LDC_bound_psib'; 

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

% tolerance on residual used in spgl1 solver
sigma = .001;

pc_solver = 1;  %0 is LS, 1 is mmv and 2 is individual spg

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

% everything from here is normalized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pc_val = 0; 

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
u_low_int = interp1(x_l, B, x_h);  
err_low = norm(u_low_int-A)/norm(A); 

A1 = psi_low(:,2:end); 
% A1 = psi_ref; 
mu = max(sum(abs(A1).^2,2)); 
%
M = 1/size(A1,1)*(A1'*A1); 
s = norm(M - eye(length(M)));
1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over N_hi, or n to compute epsilon tau, or r for reduced basis... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I want a function that does these things. 

n_vec = 35:5:100; 
% n_vec = 30:1:35; 

r = 3; 
R = r+10; 

%%% Check that q > 0 
% works for 33 for LDC - Bound is very loose at this point. 
% Doesn't seem useful in general... 

% n_vec = 3:4; 


n_samps = length(n_vec); 

err_bi_vec = zeros(1,n_samps); 
err_bound_vec = zeros(1,n_samps); 
err_ahat_vec = zeros(1,n_samps); 
err_bound_N_vec = zeros(1,n_samps); 
err_bound_n_vec = zeros(1,n_samps); 

err_T_vec = zeros(1,n_samps); 
err_E_vec = zeros(1,n_samps); 

omega_vec = zeros(1,n_samps); 
q_N_vec = zeros(1,n_samps); 
q_n_vec = zeros(1,n_samps); 

err_E_exact_vec = zeros(1,n_samps); 
err_T_exact_vec = zeros(1,n_samps); 
err_det_vec = zeros(1,n_samps); 
err_samp_vec = zeros(1,n_samps); 
err_CS_1_vec = zeros(1,n_samps); 
err_CS_2_vec = zeros(1,n_samps); 


for i_samp = 1:length(n_vec)

[err_bi_vec(i_samp), err_bound_vec(i_samp), err_ahat_vec(i_samp), ...
    err_bound_N_vec(i_samp), err_bound_n_vec(i_samp), err_Bhat, ...
    err_T_vec(i_samp), err_E_vec(i_samp), omega, q_N, q_N_hi, Ir, ...
    err_E_exact_vec(i_samp),err_T_exact_vec(i_samp),err_CS_1_vec(i_samp),...
    err_CS_2_vec(i_samp)] = ...
    br_bi_bound_test(B, A, n_vec(i_samp), R, psi_ref, psi_low, c_low, sigma, ...
    r, p, xi_low, pc_solver); 
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note - using N_hi+10 samples to compute epsilon tau component of bound. 
% I should consider changing this. 

%%% Bi and bound
figure
p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p1,p2],{'Bi','Bound'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)

if save_on == 1
    saveas(gcf,strcat('bound_plots/', res_name, '_performance'),'png')
end    

%%% Bound components
% N_vec is very close to 0. -added to n
% err_ahat_vec is almost constant. 
figure
p1 = semilogy(n_vec,err_ahat_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
% p3 = semilogy(n_vec,err_bound_N_vec,'-d','Color',c4,'LineWidth',LW,'MarkerSize',MS);
p4 = semilogy(n_vec,err_bound_n_vec+err_bound_N_vec,'-v','Color',c4,'LineWidth',LW,'MarkerSize',MS);

axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
% legend([p2,p1, p3, p4],{'Total', 'Det','N','n'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
legend([p2,p1, p4],{'Total', 'Det','Samp'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('Bound Contributions','interpreter', 'latex', 'fontsize', FS_leg)

if save_on == 1
    saveas(gcf,strcat('bound_plots/', res_name, '_contributions'),'png')
end   

figure
p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(n_vec,err_ahat_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p4 = semilogy(n_vec,err_bound_n_vec+err_bound_N_vec,'-v','Color',c4,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p1,p2, p3, p4],{'Bi','Bound', 'Det', 'Samp'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)

if save_on == 1
    saveas(gcf,strcat('bound_plots/', res_name, '_all'),'png')
end    


figure
p1 = semilogy(n_vec,err_T_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_T_exact_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert \mathbf{T} \Vert$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p1,p2],{'$\epsilon(\tau)$', 'T'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('$\mathbf{T}$ Test','interpreter', 'latex', 'fontsize', FS_leg)

if save_on == 1
    saveas(gcf,strcat('bound_plots/', res_name, '_T'),'png')
end   

figure
p1 = semilogy(n_vec,err_E_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_E_exact_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert \mathbf{E} \Vert$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p1,p2],{'$\epsilon(\tau)$', 'E'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('$\mathbf{E}$ Test','interpreter', 'latex', 'fontsize', FS_leg)

if save_on == 1
    saveas(gcf,strcat('bound_plots/', res_name, '_E'),'png')
end   

CS = [err_CS_1_vec(1), err_CS_2_vec(1)]

% Plot deterministic error with standard, with exact E and T and with exact
% E and T and pre cauchy schwarz. 

1; 

% err_ahat_vec is standard. If I use exact E and T. 
error_exact_ET = 2.*err_E_exact_vec+err_T_exact_vec*err_Bhat; 
error_pre_CS = err_E_exact_vec+err_T_exact_vec*err_Bhat+err_CS_1_vec; 

figure
% p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
% % p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(n_vec,err_ahat_vec,'-d','Color',c3,'LineWidth',LW,'MarkerSize',MS);
hold on
p4 = semilogy(n_vec,error_exact_ET,'-s','Color',c4,'LineWidth',LW,'MarkerSize',MS);
p5 = semilogy(n_vec,error_pre_CS,'-^','Color',c5,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p3, p4, p5],{'$\epsilon(\tau)$', 'E, T exact', 'pre CS'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)

if save_on == 1
    saveas(gcf,strcat('bound_plots/', res_name, '_tightness'),'png')
end    
