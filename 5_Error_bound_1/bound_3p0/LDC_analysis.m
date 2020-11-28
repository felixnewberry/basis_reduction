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

% LDC 

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


% append results with new data: 

load('u_hi.mat')
u_high = u_matrix_0; 
load('u_nu_vec_hi.mat')
xi_hi_2 = xi_hi*2-1; 

u_ref = [u_ref; u_high];
xi_ref = [xi_ref; xi_hi_2]; 

n_est = 200;
n_points = length(x_h); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 4;                          % PCE order
d = 2;                          % Stochastic dimension

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
t = 0.95; 

n_reps = 30; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Uc = u_low'; Uf = u_ref'; 
B = Uc/norm(Uc,'fro');
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
%%% Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = 3; 
R = r+10;

N_hi = 20; 

err_bi_sum_rep  = zeros(n_est,n_reps);
err_bi_mean_rep = zeros(n_points, n_reps); % mean(err_bi_mat,2)
mu_rep          = zeros(1, n_reps); 
rho_k_rep       = zeros(1, n_reps); 
Y_Nh_rep        = zeros(n_points, n_reps); 
theta_rep       = zeros(n_points, n_reps); 
U_bar_rep       = zeros(n_points, n_reps); 
U_hat_rep       = zeros(n_points, n_reps); 
Y_num_rep       = zeros(n_points, n_reps); % vecnorm(Y_num, 2, 2)
Y_den_rep       = zeros(n_points, n_reps); 
err_mat_rep     = zeros(3, n_reps); 
p_39_rep        = zeros(n_points, n_reps); 
bound_40_rep    = zeros(n_points, n_reps); 
p_41_rep        = zeros(1, n_reps); 
bound_42_rep    = zeros(1, n_reps); 


for i_rep = 1:n_reps
    
% [err_bi_sum, err_bi_mat, mu,rho_k, Y_Nh_vec, theta_vec, U_bar_vec,...
%     U_hat_vec,  Y_num, Y_den, err_mat, p_39, bound_40, p_41, bound_42] = ...
%     br_ep_tau_3p0(B, A_inf, N_hi, R, psi_ref, c_low, sigma, r, n_est, t, C);

[err_bi_sum_rep(:,i_rep), err_bi_mean_rep(:,i_rep), mu_rep(i_rep), ...
    rho_k_rep(i_rep), Y_Nh_rep(:,i_rep), theta_rep(:,i_rep), ...
    U_bar_rep(:,i_rep), U_hat_rep(:,i_rep),  Y_num_rep(:,i_rep), ...
    Y_den_rep(:,i_rep), err_mat_rep(:,i_rep), p_39_rep(:,i_rep), ...
    bound_40_rep(:,i_rep), p_41_rep(i_rep), bound_42_rep(i_rep)] = ...
    br_ep_tau_3p0(B, A_inf, N_hi, R, psi_ref, c_low, sigma, r, n_est, t, C);
end


err_bi_sum = mean(err_bi_sum_rep,2); 
err_bi_mat = mean(err_bi_mean_rep,2); 
mu = mean(mu_rep); 
rho_k = mean(rho_k_rep); 
Y_Nh_vec = mean(Y_Nh_rep,2); 
theta_vec = mean(theta_rep,2); 
U_bar_vec = mean(U_bar_rep,2); 
U_hat_vec = mean(U_hat_rep,2); 
Y_num = mean(Y_num_rep,2); 
Y_den = mean(Y_den_rep,2); 
err_mat = mean(err_mat_rep,2); 
p_39 = mean(p_39_rep,2); 
bound_40 = mean(bound_40_rep,2); 
p_41 = mean(p_41_rep); 
bound_42 = mean(bound_42_rep); 

term_1 = (1+4*mu/N_hi);

% could compute without max values too... 
term_25 = theta_vec./N_hi.*U_hat_vec.^2;
term_27 = Y_Nh_vec.*theta_vec./N_hi.*U_bar_vec.^2;
term_29 = r*max(Y_Nh_vec)*max(theta_vec)./N_hi*rho_k^2;

% Y_Nh_vec2 max is quite large - 5e5... at index 50. 

bound_25 = term_1*term_25; % presently matches 27 by design
bound_27 = term_1*term_27;
bound_27_sum = term_1*sum(term_27);
bound_29 = term_1*term_29;

% something is wrong
efficacy_27 = sqrt(bound_27./mean(err_bi_mat,2));
efficacy_27_sum = sqrt(bound_27_sum./mean(err_bi_sum));
efficacy_29 = sqrt(bound_29./mean(err_bi_sum));

efficacy_40 = sqrt(bound_40./mean(err_bi_mat,2));
efficacy_42 = sqrt(bound_42./mean(err_bi_sum));


[efficacy_29, efficacy_27_sum, efficacy_42, mean(p_39), p_41] % interesting to compare... 


figure
p1 = plot(x_l, efficacy_27,'-x','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(x_l, efficacy_40,'--o','Color',c2,...
    'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Efficacy', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2],{'Efficacy 27', 'Efficacy 40'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')

if save_on == 1
    saveas(gcf,'plots/LDC_eq_27_efficacy','epsc')
    saveas(gcf,'plots/LDC_eq_27_efficacy','fig')
end

% share this
figure
plot(x_l,theta_vec,'-x','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Theta$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on == 1
    saveas(gcf,'plots/GT_theta','epsc')
    saveas(gcf,'plots/GT_theta','fig')
end

% share this
figure
semilogy(x_l, Y_Nh_vec,'-x','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$Y$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on == 1
    saveas(gcf,'plots/GT_Y','epsc')
    saveas(gcf,'plots/GT_Y','fig')
end

max(Y_Nh_vec) % should be between 0 and 1... 

% check rho_k difference: compare term from (29), (27) and (25)
% share this 
[r*rho_k^2, sum(U_hat_vec.^2./Y_Nh_vec)] 

% rho_k doesn't properly bound interpolative decomposition as expected
% (This could be because the number of sampls is small)
% It does correctly bound when I increase N_H... 
% Y remains consistently much larger than (0,1]... this is a major problem.

figure
p1 = semilogy(x_l, sqrt(mean(err_bi_mat,2)),'-o','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(x_l, sqrt(bound_27),'--s','Color',c3,...
    'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(x_l, sqrt(bound_40),'-.x','Color',c4,...
    'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% first and last points are zero so don't plot these: 
xlim([x_l(2),x_l(64)])
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2,p3],{'True Mean','Bound 27', 'Bound 40'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')
if save_on == 1
    saveas(gcf,'plots/LDC_eq_27_bound','epsc')
    saveas(gcf,'plots/LDC_eq_27_bound','fig')
end
% is there a method to have approximation improve for specific points?
% maybe this isn't crazy? - Yvec varies from 0.1 to 2e3. But the mean is
% 60, and sqrt(60) is 7.75. 
% There's really just a few points that are large ratio of 35-40
% index 28, 47 (less so) and 65 (major)
% 28 has no significance, 47 is second turning point, 65 is end... takes
% ratio 1e-46 and 1e-49... error is negligable. 
% create write up - describe this phenomenon - do analysis and N_H changes
% - and r? 
% Do airfoil. 
% Consult Alireza and Jared. LOTS to do!!! maybe just N_H but not yet r... 

% % figure 
% plot(a1)
% plot(a2) is revealing. 
% Send all this stuff to output and examine... 
% also plot the non sqrt ratio...


% Y_num, Y_den, err_mat

% err_mat is err_mat = [err_Abar, err_Ahat, rho_k]; 
err_mat

figure
p1 = plot(x_l, Y_num,'-o','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(x_l, Y_den,'-x','Color',c2,...
    'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% first and last points are zero so don't plot these: 
% xlim([2,64])
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2],{'$\Vert U_H(:,i) - \hat{U}_H(:,i) \Vert_2$',...
    '$\Vert U_H(:,i) - \bar{U}_H(:,i) \Vert_2$'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,'plots/LDC_Y_ratio','epsc')
    saveas(gcf,'plots/LDC_Y_ratio','fig')
end



              