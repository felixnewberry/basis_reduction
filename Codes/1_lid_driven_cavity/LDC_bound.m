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
% t = 0.95; 
% phi_t = normcdf(t); p_t is hmm... 
%%% a check on t influence. p(t) = .6 and .95. % 30 reps? 
% % t = 0.95 yields probabilities of 0.709... 
% t = 4; %is .877
% % t = 0.6 is .60
% t = 0.6; 
t=0.95;

n_reps = 2;  

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Efficacy loop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Where is the parfor stuff? Maybe I didn't do that for these... which is
% % silly. Try run the LDC now and see how it goes? 
% % N_hi_vec = 3:50; 
% % r_vec = 3:20; 
% 
% N_hi_vec = 15;
% r_vec = 3; 
% 
% % N_hi_vec = [3,20];
% % r_vec = [3, 50]; 
% % n_reps = 2; 
% 
% 
% length_n = length(N_hi_vec);
% length_r = length(r_vec); 
% 
% n_r_results{n_reps} = [];
% 
% % Repeat
% % parfor i_rep = 1:n_reps
% for i_rep = 1:n_reps
% 
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
%     save('Results/LDC_efficacy', 'r_vec', 'efficacy_mat', 'prob_mat', 'N_hi_vec')
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Single point:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How to speed this up? 
% Really we just want a few things out compute those in the for loop then
% send out so that this is tidier. 
r = 3; 
R = r+10;

N_hi = 15; 

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
p_33_rep        = zeros(n_points, n_reps); 
bound_34_rep    = zeros(n_points, n_reps); 
p_35_rep        = zeros(1, n_reps); 
bound_36_rep    = zeros(1, n_reps); 
U_bar_fro_rep    = zeros(1, n_reps); 
U_bar_2_rep    = zeros(1, n_reps); 


for i_rep = 1:n_reps
   

[err_bi_sum_rep(:,i_rep), err_bi_mean_rep(:,i_rep), mu_rep(i_rep), ...
    rho_k_rep(i_rep), Y_Nh_rep(:,i_rep), theta_rep(:,i_rep), ...
    U_bar_rep(:,i_rep), U_hat_rep(:,i_rep),  Y_num_rep(:,i_rep), ...
    Y_den_rep(:,i_rep), err_mat_rep(:,i_rep), p_33_rep(:,i_rep), ...
    bound_34_rep(:,i_rep), p_35_rep(i_rep), bound_36_rep(i_rep), ...
    U_bar_fro_rep(i_rep), U_bar_2_rep(i_rep)] = ...
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
p_33 = mean(p_33_rep,2); 
bound_34 = mean(bound_34_rep,2); 
p_35 = mean(p_35_rep); 
bound_36 = mean(bound_36_rep); 
U_bar_fro = mean(U_bar_fro_rep); 
U_bar_2 = mean(U_bar_2_rep); 

term_1 = (1+4*mu/N_hi);

% 25 to 23. 
% 27 to A6. 
% 39 to 33
% 41 to 35
% 40 to 34
% 42 to 36

% could compute without max values too... 
term_23 = theta_vec./N_hi.*U_hat_vec.^2;
term_A6 = Y_Nh_vec.*theta_vec./N_hi.*U_bar_vec.^2;
term_29 = r*max(Y_Nh_vec)*max(theta_vec)./N_hi*rho_k^2;

% Y_Nh_vec2 max is quite large - 5e5... at index 50. 

bound_23 = term_1*term_23; 
bound_A6 = term_1*term_A6;
bound_A6_sum = term_1*sum(term_A6);
bound_A8 = term_1*term_29;

% something is wrong
efficacy_23 = sqrt(bound_23./mean(err_bi_mat,2));
efficacy_A6 = sqrt(bound_A6./mean(err_bi_mat,2));
efficacy_A6_sum = sqrt(bound_A6_sum./mean(err_bi_sum));
efficacy_A8 = sqrt(bound_A8./mean(err_bi_sum));

efficacy_34 = sqrt(bound_34./mean(err_bi_mat,2));
efficacy_36 = sqrt(bound_36./mean(err_bi_sum));


efficacy_vec = [efficacy_A8, efficacy_A6_sum, efficacy_36, mean(p_33), p_35];
% interesting to compare... 

% check rho_k difference: compare term from (29?), (A6) and (23)
% share this 

% U_bar_fro^2 = sum(U_bar_vec.^2)

% rho_vec = [r*rho_k^2, sum(U_hat_vec.^2./Y_Nh_vec)]; 
rho_vec = [sqrt(r*rho_k^2), sqrt(sum(U_bar_vec.^2))]; 
rho_vec_2 = [sqrt(r*rho_k^2), U_bar_2]; 
% rho_vec = [sqrt(r*rho_k^2), sqrt(sum(U_hat_vec.^2./Y_Nh_vec))]; 
% rho_vec = [sqrt(r*rho_k^2), U_bar_fro]; 
% Ah, I see, I've undone the change from hat to bar... this may not be
% correct... 


% Is the second term here correct?? 
% rho_k doesn't properly bound interpolative decomposition as expected
% (This could be because the number of sampls is small)
% It does correctly bound when I increase N_H... 
% Y remains consistently much larger than (0,1]... this is a major problem.

% if save_on == 1
%     save('Results/LDC_bound_results', 'r_vec', 'efficacy_mat',...
%     'N_hi_vec', 'r', 'R', 'N_hi', 'n_reps', 'x_h', ...
%     'efficacy_27', 'efficacy_40', 'theta_vec', 'Y_Nh_vec', 'err_bi_mat',...
%     'bound_27', 'bound_40', 'efficacy_vec', 'rho_vec');
% end
% if save_on == 1
%     save('Results/LDC_bound_results_test_t_p6_p60', 'r_vec', 'efficacy_mat',...
%     'N_hi_vec', 'r', 'R', 'N_hi', 'n_reps', 'x_h', ...
%     'efficacy_23', 'efficacy_A6', 'efficacy_34', 'theta_vec', 'Y_Nh_vec', 'err_bi_mat',...
%     'bound_23','bound_A6', 'bound_34', 'efficacy_vec', 'rho_vec');
% end
if save_on == 1
    save('Results/LDC_bound_results_test_t', 'r_vec', 'efficacy_mat',...
    'N_hi_vec', 'r', 'R', 'N_hi', 'n_reps', 'x_h', ...
    'efficacy_23', 'efficacy_A6', 'efficacy_34', 'theta_vec', 'Y_Nh_vec', 'err_bi_mat',...
    'bound_23','bound_A6', 'bound_34', 'efficacy_vec', 'rho_vec');
end

% check p39 and p41
