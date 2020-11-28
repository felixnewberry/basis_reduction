clear all
close all
clc

%%% Sampling bound check. 

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

res_name = 'Airfoil_bound'; 

load lowFiResults_6d2.mat
load highFiResults_6d2.mat;

load Inputs_6d2.mat;
xi_ref =[(1/3)*alpha',(4*c-4)',(50*m-1)',(10*(mach-.2)-2)',(10*p -4)',((t - .1)*50 - 1)'];

xi_low = xi_ref; 

% Grid points of interest
gridpt_l = 1:size(lowFiResults(:,1)); 
gridpt_h = 1:size(highFiResults(:,1)); 

% u should be number of samples x spatial/temporal dimension of the
% problem. size: n_samples x n_gridpoints
u_ref = highFiResults(gridpt_h,:)';
u_low = lowFiResults(gridpt_l,:)';


% What if I aproximate just 200 samples total? Will this change anything? 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 5;                          % PCE order
d = 6;                          % Stochastic dimension

% tolerance on residual used in spgl1 solver
sigma = .001;

pc_solver = 0;  %0 is LS, 1 is mmv and 2 is individual spg

n_sim = length(xi_ref); 

t_setup = toc; 
fprintf('Setup with time : %d s.\n', t_setup);
timerval = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Max 1200. 
data_samples = 1:1200; 

xi_low = xi_low(data_samples,:); 
xi_ref = xi_ref(data_samples,:); 

Uc = u_ref(data_samples,:)'; 
% Uc = u_low(data_samples,:)'; 

Uf = u_ref(data_samples,:)'; 

B = Uc/norm(Uc,'fro');
A = Uf/norm(Uf,'fro');

% everything from here is normalized
rng(42); 
% samples are repeatably random. 


load('err_bi_Uf_10_samps'); 
n_vec = [50,100,200, 400, 600, 800]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pc_val = 1; 
    
%%% reference solution
[c_ref, psi_ref, err_val] = my_pce(xi_ref, p, A', sigma, pc_solver, pc_val); 
t_ref = toc(timerval); 
fprintf('Reference solution : %d s.\n', t_ref);
timerval = tic; 

%%% Low-fidelity solution
[c_low, psi_low, err_val_low] = my_pce(xi_low, p, B', sigma, pc_solver, pc_val); 
t_low = toc(timerval); 
% t1 = toc
% t_ref
fprintf('Low fidelity solution : %d s.\n', t_low);

%%% what am I doing wrong here? 
% A1 = psi_ref(:,2:end); 
A1 = psi_ref/sqrt(n_sim); 
mu = max(sum(abs(A1).^2,2)); 
%
M = (A1'*A1); 
s = norm(M - eye(length(M)))
% > 1/2 - does not satisfy condition coherence motivated sampling. 
%

% Lemma 1: 
size(A1,1)*sqrt(size(A1,2)); %>= 20
Prob = 2*size(A1,2)*exp(-0.1*size(A1,1)*mu^-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Low-fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % n_vec = 200; 
% 
% % n_vec = [50,100,200, 400, 600, 800]; 
% n_vec = [20:10:100]; 
% 
% r = 4; 
% % r = 8; 
% R = r+10; % 
% % R = 20; 
% 
% n_reps = 1; 
% 
% n_samps = length(n_vec); 
% 
% err_bi_vec = zeros(n_reps,n_samps); 
% err_bound_vec = zeros(n_reps,n_samps); 
% err_ahat_vec = zeros(n_reps,n_samps); 
% err_bound_N_vec = zeros(n_reps,n_samps); 
% err_bound_n_vec = zeros(n_reps,n_samps); 
% 
% err_T_vec = zeros(n_reps,n_samps); 
% err_E_vec = zeros(n_reps,n_samps); 
% 
% omega_vec = zeros(n_reps,n_samps); 
% q_N_vec = zeros(n_reps,n_samps); 
% q_n_vec = zeros(n_reps,n_samps); 
% 
% err_E_exact_vec = zeros(n_reps,n_samps); 
% err_T_exact_vec = zeros(n_reps,n_samps); 
% err_det_vec = zeros(n_reps,n_samps); 
% err_samp_vec = zeros(n_reps,n_samps); 
% err_CS_1_vec = zeros(n_reps,n_samps); 
% err_CS_2_vec = zeros(n_reps,n_samps); 
% 
% err_bi_points_mean = zeros(128,n_samps); 
% err_bi_points_med = zeros(128,n_samps); 
% 
% for i_rep = 1:n_reps
%     
% for i_samp = 1:length(n_vec)
% 
% [err_bi_vec(i_rep, i_samp), err_bound_vec(i_rep, i_samp), ...
%     err_ahat_vec(i_rep, i_samp), err_bound_N_vec(i_rep, i_samp), ...
%     err_bound_n_vec(i_rep, i_samp), err_Bhat, err_T_vec(i_rep, i_samp), ...
%     err_E_vec(i_rep, i_samp), omega, q_N, q_N_hi, Ir, ...
%     err_E_exact_vec(i_rep, i_samp),err_T_exact_vec(i_rep, i_samp), ...
%     err_CS_1_vec(i_rep, i_samp), err_CS_2_vec(i_rep, i_samp), Ahat,...
%     err_bi_points_mean_temp, err_bi_points_med_temp] = ...
%     br_bi_bound_test(B, A, n_vec(i_samp), R, psi_ref, c_low, ...
%     sigma, r); 
% 
%     err_bi_points_mean(:,i_samp) = err_bi_points_mean(:,i_samp)+err_bi_points_mean_temp; 
%     err_bi_points_med(:,i_samp) = err_bi_points_med(:,i_samp)+err_bi_points_med_temp; 
% 
% end
% 
% end
% 
% figure
% plot(n_vec, err_bi_vec)
% % err_bi_vec is on the order of 5%... ie 5e-2 order 1$ 
% 
% % account for samples
% err_bi_points_mean = err_bi_points_mean/n_reps; 
% err_bi_points_med = err_bi_points_med/n_reps; 
% 
% figure
% plot(err_bi_points_med)

% % relative error point by point
% err_bi_points_med_rel = err_bi_points_med./vecnorm(A, 2, 2);
% 
% figure
% plot(A(40,:), 'rx')
% hold on
% plot(Ahat(40,:), 'bo')
% 
% figure
% plot(A(:,100), 'rx')
% hold on
% plot(Ahat(:,100), 'bo')
% 
% e_test = sqrt((A(:,100)- Ahat(:,100)).^2)./sqrt(mean(Ahat(:,100).^2));
% Looks good. 

% If I look at just one point. say 40. 


% med_test = median(abs(Ahat - A),2);
% 
% save('err_bi_points_2', 'err_ahat_vec', 'err_bi_vec', 'err_bi_points_med', 'err_bi_points_mean','n_vec', 'r', 'R')

% 
%  
% 
% save('err_bi_vec_3', 'err_Ahat_vec', 'err_bi_vec', 'n_vec', 'r', 'R', 'err_ahat_vec')

% figure
% plot(n_vec, err_bi_points_med, 'r')
% hold on
% plot(n_vec, err_bi_points_mean, 'b')
% order 0 to 2.5 e-5 
1; % check first plot of error bound - then show process - identify slope r - go to points. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 10 repetitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% load('err_bi_vec_1.mat')
% load('err_bi_vec_2.mat')
% load('err_bi_vec_3.mat')
% 
% figure
% p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
% hold on
% % p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
% p3 = semilogy(n_vec,err_ahat_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
% % p4 = semilogy(n_vec,err_bound_n_vec+err_bound_N_vec,'-v','Color',c4,'LineWidth',LW,'MarkerSize',MS);
% hold off
% axis tight
% xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% % grid on
% set(gcf,'Position',size_1)
% % legend([p1,p2, p3, p4],{'Bi','Bound', 'Det', 'Samp'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
% % legend([p1(1),p2, p3(1)],{'Bi','Total Bound', 'Det Bound',},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optimal r/q value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('err_bi_points_2.mat')

err_bi_mean = mean(err_bi_vec); 

err_bi_points_med_mean = mean(err_bi_points_med,1);
err_bi_points_mean_mean = mean(err_bi_points_mean,1);

err_bi_points_med2 = err_bi_points_med(1,:);
err_bi_points_mean2 = err_bi_points_mean(1,:);

%%% demonstrate appropriate error to investigate. 
figure
p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_ahat_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error $\Vert H - \hat{H} \Vert$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p1(1),p2(1)],{'Bi','Bound Det'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')


figure
p1 = plot(log(n_vec), log(err_bi_points_med.^2),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(log(n_vec), log(err_bi_points_mean.^2),'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$\log{n}$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\log{\Vert e \Vert^2}$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Optimal r calculation','interpreter', 'latex', 'fontsize', FS_leg)
legend([p1(1),p2(1)],{'mean','med'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

figure
p1 = plot(log(n_vec), log(err_bi_points_med2.^2),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(log(n_vec), log(err_bi_points_mean2.^2),'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$\log{n}$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\log{\Vert e \Vert^2}$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Optimal r calculation','interpreter', 'latex', 'fontsize', FS_leg)
legend([p1(1),p2(1)],{'mean(mean)','mean(med)'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')




%%% Calculate r for each point. 



% Try mean, median and see if same as taking average first. 

r_opt_vec = zeros(size(err_bi_points_med,1), 2); 
for i_points = 1:size(err_bi_points_med,1)
    coefs = polyfit(log(n_vec), log(err_bi_points_med(i_points,:).^2), 1);
    r_opt(i_points,1) = -coefs(1);
    coefs = polyfit(log(n_vec), log(err_bi_points_mean(i_points,:).^2), 1);
    r_opt(i_points,2) = -coefs(1);     
end

coefs = polyfit(log(n_vec), log(err_bi_points_med_mean.^2), 1);
r_med = -coefs(1);
coefs = polyfit(log(n_vec), log(err_bi_points_mean_mean.^2), 1);
r_mean = -coefs(1);

n_hist = 20; 

figure
hold on
h1 = histogram(r_opt(:,1),n_hist,'FaceColor',c1);
h2 = histogram(r_opt(:,2),n_hist,'FaceColor',c2);
hold off
legend([h1,h2],{'Median','Mean'},'interpreter', 'latex', 'fontsize', FS_leg)
xlabel('$q$ $(r)$','interpreter','latex','Fontsize',FS)
ylabel('Frequency','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
set(gcf,'Position',size_1)

% max and minimum r
r_max = max(r_opt(:)); 
r_min = min(r_opt(:)); 


% Now - testing the bound performance with r_min and r_max. 

% 1.6 in paper. 

% I can calculate L from data - assume access to all Uf - in practice just
% n samples.  - Use A rather than Uf. 
% We know n
% Calculate epsilon from r
% How to estimate the best approximation error? 

L_vec = max(abs(A),[],2); 
% q = r_max; 


q = r_max; % 0.1960


kappa = (3*log(3/2)-1)/(2+2*q); 
epsilon = 4*kappa./log(n_vec); 

% best approximation error, take minimum? 
% e_mf2 = min(err_bi_points_med(:,end)); 
% Should really be different for each point. 
e_mf2 = err_bi_points_med(:,end).^2; 
% This best approximation error, once squared, is very small. Term 1 is
% small - but same order as actual error squared... 
% The error using the best r is very loose. 

% 128 points - 9 different n
term_1 = (1+epsilon).*e_mf2; 
term_2 = 8*L_vec.^2*n_vec.^(-q); 
bound = term_1 + term_2; 

% bound barely changes with n. Therfore just plot the tightest bound for
% all points. 


figure
p1 = plot(bound(:,end),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(err_bi_points_med(:,end).^2,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = plot(term_1(:,end),'-.s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p4 = plot(term_2(:,end),'-.d','Color',c4,'LineWidth',LW,'MarkerSize',MS);

axis tight
xlabel('point index', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert e \Vert^2$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1(1),p2(1), p3(1), p4(1)],{'Point Bound','Point Error', 'Term 1', 'Term 2'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

% e_test = mean(abs(A - Ahat),2); 
% plot(e_test.^2)
% title('Optimal r calculation','interpreter', 'latex', 'fontsize', FS_leg)


%%% Tidy up results and make plots nice. Check I didn't miss something
%%% important. 

% The real errors are still very small. 
% Can I adjust this problem so that the relative errors are order 5%, and
% the errors squared are 0.25%? 

% Try do that. 


%%% r robust test: 


L_vec = max(abs(A),[],2); 
e_mf2 = err_bi_points_med(:,end).^2; 


% loop over candidate r. Check bound is satisifed... - by minimum amount 
% start with r = 0.01 and multiply by 2? 

% how to converge best? 
q_cand = [0.01, 0.1, 0.5,1, 1.5, 2, 2.5, 2.6, 2.75, 3, 3.5]; 
efficacy_med_vec = zeros(1,length(q_cand)); 

for i_q = 1:length(q_cand)
    q = q_cand(i_q); 
    
    kappa = (3*log(3/2)-1)/(2+2*q); 
    epsilon = 4*kappa./log(n_vec); 
    
    term_1 = (1+epsilon).*e_mf2; 
    term_2 = 8*L_vec.^2*n_vec.^(-q); 
    bound = term_1 + term_2; 
    
    % Upper bound?
    efficacy_mat = bound./(err_bi_points_med.^2);
    if min(efficacy_mat(:)) <= 1
        efficacy_med_vec(i_q) = NaN; 
    else
        efficacy_med_vec(i_q) = median(efficacy_mat(:)); 
    end
end


[best_med, i_q_best] = min(efficacy_med_vec);

q_best = q_cand(i_q_best)
kappa = (3*log(3/2)-1)/(2+2*q_best); 
epsilon = 4*kappa./log(n_vec); 

term_1 = (1+epsilon).*e_mf2; 
term_2 = 8*L_vec.^2*n_vec.^(-q_best); 
bound = term_1 + term_2; 

figure

subplot(1,2,1)
p1 = plot(bound(:,1),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(err_bi_points_med(:,1).^2,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = plot(term_1(:,1),'-.s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p4 = plot(term_2(:,1),'-.d','Color',c4,'LineWidth',LW,'MarkerSize',MS);

axis tight
xlabel('point index', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert e \Vert^2$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1(1),p2(1), p3(1), p4(1)],{'Point Bound','Point Error', 'Term 1', 'Term 2'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('$n = 20$ samples','interpreter', 'latex', 'fontsize', FS_leg)

subplot(1,2,2)
p1 = plot(bound(:,end),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(err_bi_points_med(:,end).^2,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = plot(term_1(:,end),'-.s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p4 = plot(term_2(:,end),'-.d','Color',c4,'LineWidth',LW,'MarkerSize',MS);

axis tight
xlabel('point index', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Vert e \Vert^2$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1(1),p2(1), p3(1), p4(1)],{'Point Bound','Point Error', 'Term 1', 'Term 2'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('$n = 100$ samples','interpreter', 'latex', 'fontsize', FS_leg)

set(gcf, 'Position', size_2)

  