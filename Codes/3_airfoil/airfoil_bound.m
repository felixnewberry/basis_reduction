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
t = 0.95; 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Efficacy loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where is the parfor stuff? Maybe I didn't do that for these... which is
% silly. Try run the LDC now and see how it goes? 
% N_hi_vec = 3:50; 
% r_vec = 3:20; 

N_hi_vec = 30; 
r_vec = 8; 

%r_vec = [3, 20]; 
%N_hi_vec = [3,50];
%n_reps = 2;

length_n = length(N_hi_vec);
length_r = length(r_vec); 

n_r_results{n_reps} = [];

% Repeat
%parfor i_rep = 1:n_reps
for i_rep = 1:n_reps
    i_rep
    n_r_results{i_rep}.efficacy = zeros(length_r,length_n); 
    n_r_results{i_rep}.prob = zeros(length_r,length_n); 

    for i_n = 1:length_n
        for i_r = 1:length_r

        [n_r_results{i_rep}.efficacy(i_r, i_n), n_r_results{i_rep}.prob(i_r, i_n)] = ...
    br_bound_practical(A_inf, N_hi_vec(i_n), psi_ref, c_low, sigma, r_vec(i_r), n_est, t, C);

        end
    end
end
1;

% Calculate mean over n_reps
stat_struct = cat(3,n_r_results{:});
efficacy_mat = mean(cat(3,stat_struct.efficacy),3);
prob_mat = mean(cat(3,stat_struct.prob),3);

% efficacy_mat = mean_ep_tau_bound./mean_bi_err; 

if save_on == 1
    save('Results/Airfoil_efficacy', 'r_vec', 'efficacy_mat', 'prob_mat', 'N_hi_vec')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single point:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = 8; 
R = r+10;

N_hi = 30; 

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
    br_bound_complete(B, A_inf, N_hi, R, psi_ref, c_low, sigma, r, n_est, t, C);
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

efficacy_vec = [efficacy_29, efficacy_27_sum, efficacy_42, mean(p_39), p_41];
% interesting to compare... 

% check rho_k difference: compare term from (29), (27) and (25)
% share this 
rho_vec = [r*rho_k^2, sum(U_hat_vec.^2./Y_Nh_vec)]; 

if save_on == 1
save('Results/Airfoil_bound_results', 'r_vec', 'efficacy_mat',...
    'N_hi_vec', 'r', 'R', 'N_hi', 'n_reps', 'x_l', ...
    'efficacy_27', 'efficacy_40', 'theta_vec', 'Y_Nh_vec', 'err_bi_mat',...
    'bound_27', 'bound_40', 'efficacy_vec', 'rho_vec');
end

% figure
% plot(x_l, efficacy_27,'-x','Color',c1,...
%     'LineWidth',LW,'MarkerSize',MS);
% axis tight
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Efficacy Equantion (27)', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% % ylim([1e-8,1])
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% if save_on == 1
%     saveas(gcf,'plots/Airfoil_eq_27_efficacy','epsc')
%     saveas(gcf,'plots/Airfoil_eq_27_efficacy','fig')
% end
% 
% n_hist = 20; 
% 
% % No need to share
% figure
% h1 = histogram(efficacy_27,n_hist,'FaceColor',c1);
% xlabel('$efficacy_{27}$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% % axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% % No need to share
% figure
% h1 = histogram(theta_vec,n_hist,'FaceColor',c1);
% xlabel('$\Theta$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% % axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% % share this
% figure
% plot(x_l,theta_vec,'-x','Color',c1,...
%     'LineWidth',LW,'MarkerSize',MS);
% axis tight
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$\Theta$', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% % ylim([1e-8,1])
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% if save_on == 1
%     saveas(gcf,'plots/Airfoil_theta','epsc')
%     saveas(gcf,'plots/Airfoil_theta','fig')
% end
% 
% % share this
% figure
% plot(x_l, Y_Nh_vec,'-x','Color',c1,...
%     'LineWidth',LW,'MarkerSize',MS);
% axis tight
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$Y$', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% % ylim([1e-8,1])
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% if save_on == 1
%     saveas(gcf,'plots/Airfoil_Y','epsc')
%     saveas(gcf,'plots/Airfoil_Y','fig')
% end
% 
% % no need to share this
% figure
% h1 = histogram(log(Y_Nh_vec), n_hist,'FaceColor',c1);
% xlabel('$log(Y)$','interpreter','latex','Fontsize',FS)
% ylabel('Frequency','interpreter','latex','Fontsize',FS)
% % axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% max(Y_Nh_vec) % should be between 0 and 1... 
% 
% % check rho_k difference: compare term from (29), (27) and (25)
% % share this 
% [r*rho_k^2, sum(U_hat_vec.^2./Y_Nh_vec)] 
% 
% % rho_k doesn't properly bound interpolative decomposition as expected
% % (This could be because the number of sampls is small)
% % It does correctly bound when I increase N_H... 
% % Y remains consistently much larger than (0,1]... this is a major problem.
% 
% figure
% p1 = semilogy(x_l, sqrt(mean(err_bi_mat,2)),'-o','Color',c1,...
%     'LineWidth',LW,'MarkerSize',MS);
% hold on
% p2 = semilogy(x_l, sqrt(bound_27),'--s','Color',c3,...
%     'LineWidth',LW,'MarkerSize',MS);
% p3 = semilogy(x_l, sqrt(bound_40),'-.x','Color',c4,...
%     'LineWidth',LW,'MarkerSize',MS);
% hold off
% axis tight
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% % first and last points are zero so don't plot these: 
% % xlim([x_l(2),x_l(64)])
% % ylim([1e-8,1])
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% legend([p1,p2,p3],{'True Mean','Bound 27', 'Bound 40'}...
%     ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')
% 
% if save_on == 1
%     saveas(gcf,'plots/Airfoil_eq_27_bound','epsc')
%     saveas(gcf,'plots/Airfoil_eq_27_bound','fig')
% end
% 
% % is there a method to have approximation improve for specific points?
% % maybe this isn't crazy? - Yvec varies from 0.1 to 2e3. But the mean is
% % 60, and sqrt(60) is 7.75. 
% % There's really just a few points that are large ratio of 35-40
% % index 28, 47 (less so) and 65 (major)
% % 28 has no significance, 47 is second turning point, 65 is end... takes
% % ratio 1e-46 and 1e-49... error is negligable. 
% % create write up - describe this phenomenon - do analysis and N_H changes
% % - and r? 
% % Do airfoil. 
% % Consult Alireza and Jared. LOTS to do!!! maybe just N_H but not yet r... 
% 
% % % figure 
% % plot(a1)
% % plot(a2) is revealing. 
% % Send all this stuff to output and examine... 
% % also plot the non sqrt ratio...
% 
% 
% % Y_num, Y_den, err_mat
% 
% % err_mat is err_mat = [err_Abar, err_Ahat, rho_k]; 
% err_mat
% 
% figure
% p1 = plot(x_l, vecnorm(Y_num, 2, 2),'-o','Color',c1,...
%     'LineWidth',LW,'MarkerSize',MS);
% hold on
% p2 = plot(x_l,vecnorm(Y_den, 2, 2),'-x','Color',c2,...
%     'LineWidth',LW,'MarkerSize',MS);
% hold off
% axis tight
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% % first and last points are zero so don't plot these: 
% % xlim([2,64])
% % ylim([1e-8,1])
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% legend([p1,p2],{'$\Vert U_H(:,i) - \hat{U}_H(:,i) \Vert_2$',...
%     '$\Vert U_H(:,i) - \bar{U}_H(:,i) \Vert_2$'}...
%     ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
% 
% if save_on == 1
%     saveas(gcf,'plots/Airfoil_Y_ratio','epsc')
%     saveas(gcf,'plots/Airfoil_Y_ratio','fig')
% end
% 
% % shape is right - but the bound is wrong... where are things awry?
% % examine. 
% 
% % theta_vec looks decent
% % Y_Nh_vec2 also looks good - though it's not very large. Max value is
% % 0.05. 
% 
% %%% sampling less then r doesn't make much sense... what to do?
% % N_hi_vec = [5,6,7,8,9,10:10:200]; 
% % 
% % n_rep = 10; 
% % 
% % err_bi = zeros(n_rep,length(N_hi_vec));
% % err_ep_tau_bound = zeros(1,length(N_hi_vec));
% % err_Bhat = zeros(1,length(N_hi_vec));
% % err_E = zeros(1,length(N_hi_vec));
% % err_T = zeros(1,length(N_hi_vec));
% % err_E_exact = zeros(1,length(N_hi_vec));
% % err_T_exact = zeros(1,length(N_hi_vec));
% % err_CS_pre = zeros(1,length(N_hi_vec));
% % err_CS_post = zeros(1,length(N_hi_vec));
% % err_bi_N = zeros(n_rep,length(N_hi_vec));
% % err_bi_inf_n = zeros(n_rep,length(N_hi_vec));
% % err_bi_inf_N = zeros(n_rep,length(N_hi_vec));
% 
% % for i_rep =1:n_rep
% %     for i_N_hi = 1:length(N_hi_vec)
% %         [err_bi(i_rep, i_N_hi), err_ep_tau_bound(i_N_hi), err_Bhat(i_N_hi),...
% %             err_E(i_N_hi), err_T(i_N_hi), Ir, err_E_exact(i_N_hi),...
% %             err_T_exact(i_N_hi), err_CS_pre(i_N_hi), err_CS_post(i_N_hi), err_bi_N(i_rep,i_N_hi), err_bi_inf_N(i_rep,i_N_hi), err_bi_inf_n(i_rep,i_N_hi)] = ...
% %                     br_ep_tau_error_dev_Nrand(B, A_inf, N_hi_vec(i_N_hi), R, psi_ref, psi_low, c_low, sigma, ...
% %                     r, p, xi_low, pc_solver, n_est);
% %     end
% % end
% 
% % note that n samps in bound is fixed. Right now I'm changing n high
% % samples in the bi estimate. ie R is fixed n is changing. 


              
