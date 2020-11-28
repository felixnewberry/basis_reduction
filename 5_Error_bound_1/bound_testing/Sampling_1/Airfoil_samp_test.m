% Sampling Error

% Fix up 

% find psi_bi based on high with many r. 
% do this for airfoil. 

% Goal is to ensure almost all error is derived from sampling. Psi should
% be great. 
% Not currently working. 
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

%% Low_fidelity 
% u_low_int = interp1(x_l, B, x_h);  
% err_low = norm(B-A); 
% % err_val_low % p5, 0.001, p4 is 9.0e-4, p3 is 0.0044
% 1; 
% 
% %%% use all high-fidelity
% c_low = c_ref; 
% psi_low = psi_ref; 


% calibrate PC (for high-fid data)
% PC validation error
% p = 2: 2.6e-2
% p = 3: 4.5e-3
% p = 4: 4.9e-3
% p = 5: 7.9e-3 %%
% p = 6: 1.7e-2
% p = 7: 4.6e-3
% p = 8: 4.2e-3
% p = 9: 3.7e-3

% Reduce truncation error, set high r. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over N_hi, or n to compute epsilon tau, or r for reduced basis... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_vec = 50:10:100; 
% n_vec = 50:50:1000; 
% n_vec = [50,100,200, 400, 600, 800]; 
% n_vec = [500]; 
% n_vec = 500; 
% n_vec = 600; 
n_vec = 200; 

% what if n samples to solve cor c_bi exceeded the N that are estimated? 
% Then comparing n to N samples doesn't make proper sense? No, it does when
% comparing to infinity. Ponder some. 
% n_vec = 20:21; 

r = 30; 
% r = 8; 
% R = r+10; % ensure tiny truncation error - % check. 
R = r+10; % ensure tiny truncation error - % check. 
% R = 20; 

n_reps = 1; 

n_samps = length(n_vec); 

err_bi_vec = zeros(n_reps,n_samps); 
err_bound_vec = zeros(n_reps,n_samps); 
err_ahat_vec = zeros(n_reps,n_samps); 
err_bound_N_vec = zeros(n_reps,n_samps); 
err_bound_n_vec = zeros(n_reps,n_samps); 

err_T_vec = zeros(n_reps,n_samps); 
err_E_vec = zeros(n_reps,n_samps); 

omega_vec = zeros(n_reps,n_samps); 
q_N_vec = zeros(n_reps,n_samps); 
q_n_vec = zeros(n_reps,n_samps); 

err_E_exact_vec = zeros(n_reps,n_samps); 
err_T_exact_vec = zeros(n_reps,n_samps); 
err_det_vec = zeros(n_reps,n_samps); 
err_samp_vec = zeros(n_reps,n_samps); 
err_CS_1_vec = zeros(n_reps,n_samps); 
err_CS_2_vec = zeros(n_reps,n_samps); 

for i_rep = 1:n_reps
    
for i_samp = 1:length(n_vec)

[err_bi_vec(i_rep, i_samp), err_bound_vec(i_rep, i_samp), ...
    err_ahat_vec(i_rep, i_samp), err_bound_N_vec(i_rep, i_samp), ...
    err_bound_n_vec(i_rep, i_samp), err_Bhat, err_T_vec(i_rep, i_samp), ...
    err_E_vec(i_rep, i_samp), omega, q_N, q_N_hi, Ir, ...
    err_E_exact_vec(i_rep, i_samp),err_T_exact_vec(i_rep, i_samp), ...
    err_CS_1_vec(i_rep, i_samp), err_CS_2_vec(i_rep, i_samp)] = ...
    br_bi_bound_test(B, A, n_vec(i_samp), R, psi_ref, c_low, ...
    sigma, r); 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note - using N_hi+10 samples to compute epsilon tau component of bound. 
% I should consider changing this. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Curious behavior!! Figure out 
[err_bi_vec, err_bound_vec, err_ahat_vec, err_bound_n_vec+err_bound_N_vec]

% Using all 1200 samples of epsilon tau, and 500 for sampling error. 
% Bi-fidelity is larger than the bound 4.5e-3 vs 3.6e-3... 
% How is this the case?

% If I use low-fidelity... ? It appears to work. 
% Does not work for high-fidelity dataset. Something must be wrong. What is the comparison?? 
% Low-fid: bi is 0053, bound 0306.
% High-fid: bi is 0038, bound 0029.

%Hmm. 

% With 200 data samples total. 
% L: bi is 0.0035, bound is 0.0310
% H: bi and then bound
%    1.0e-03 *
% 
%     0.8164
%     0.6663

% I should test an r that is not so crazy? 

% I need to bring either the bound up, or the bi down. Don't know how. Is
% theory wrong? Examine. 

% How about I make the sampling random - do the same seed, but have it
% random. Check that results are consistent. This may be at several layers.
% Figure out. 

% Reset - plot error vs n_samples - do 10 reps. map area first. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('err_bi_Uf_10_samps'); 
n_vec = [50,100,200, 400, 600, 800]; 

figure
p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(n_vec,err_ahat_vec,'-o','Color',c3,'LineWidth',LW,'MarkerSize',MS);
% p4 = semilogy(n_vec,err_bound_n_vec+err_bound_N_vec,'-v','Color',c4,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
% legend([p1,p2, p3, p4],{'Bi','Bound', 'Det', 'Samp'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
legend([p1(1),p2, p3(1)],{'Bi','Total Bound', 'Det Bound',},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

title('Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)
ylim([2.9e-3, 1e-2])

% err_ahat_vec is standard. If I use exact E and T. 
error_exact_ET = 2.*err_E_exact_vec+err_T_exact_vec*err_Bhat; 
error_pre_CS = err_E_exact_vec+err_T_exact_vec*err_Bhat+err_CS_1_vec; 

% figure
% % p1 = semilogy(n_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
% hold on
% % % p2 = semilogy(n_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
% p3 = semilogy(n_vec,err_ahat_vec,'-d','Color',c3,'LineWidth',LW,'MarkerSize',MS);
% hold on
% p4 = semilogy(n_vec,error_exact_ET,'-s','Color',c4,'LineWidth',LW,'MarkerSize',MS);
% p5 = semilogy(n_vec,error_pre_CS,'-^','Color',c5,'LineWidth',LW,'MarkerSize',MS);
% hold off
% axis tight
% xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% % grid on
% set(gcf,'Position',size_1)
% legend([p3, p4, p5],{'$\epsilon(\tau)$', 'E, T exact', 'pre CS'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
% title('Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)
% 
% % save('err_bi_Uf_10_samps', 'err_bi_vec')