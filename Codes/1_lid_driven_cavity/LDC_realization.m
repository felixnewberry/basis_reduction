clear all
close all
clc

%%% Lid driven cavity

% Realization

tic

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

save_on = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I should re-run this with 64x64 grid as high, and try 8x8 as well as 4x4
% as low. 

results_file_name = 'LDC_results'; 

% load lowFiResults_LDC
% load xi_low_LDC
% 
% load highFiResults_LDC.mat
% load xi_high_LDC.mat
% 
% xi_ref = xi_high; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% New data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('u_64_f_2.mat')
% 200 samples. 65 points
highFiResults = u_matrix_0'; 

Uc_nom_u = load('Nom_u_mid.mat', 'Uc','Ub','sb');
lowFiResults = Uc_nom_u.Uc; 

load 'x_64.mat'
x_h = x_64(:,1); 
x_l = x_h; 

load('xi_mat_2.mat')
xi_ref = xi_2; 
xi_low = xi_ref; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Grid points of interest
gridpt_l = 1:size(lowFiResults(:,1)); 
gridpt_h = 1:size(highFiResults(:,1)); 

% u should be number of samples x spatial/temporal dimension of the
% problem. size: n_samples x n_gridpoints
u_ref = highFiResults(gridpt_h,:)';
u_low = lowFiResults(gridpt_l,:)';

% %%% Low and high coordinates
% % n_cell_l = 8;
% n_cell_l = 4;
% x_l = linspace(0,1,n_cell_l+1);
% x_l = 0.5*(cos(pi*(x_l-1)/2)+1);
% 
% % n_cell_h = 32;
% n_cell_h = 62;
% x_h = linspace(0,1,n_cell_h+1);
% x_h = 0.5*(cos(pi*(x_h-1)/2)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do n = 25, r = 3. 

p = 4;                          % PCE order
d = 2;                          % Stochastic dimension

N_hi = 15;      % Number high-fidelity samples

r = 3;                  % KL order

% tolerance on residual used in spgl1 solver
sigma = .001;

% Number of repetitions
n_reps = 1; 

pc_solver = 2;  %0 is LS, 1 is mmv and 2 is individual spg

t_setup = toc; 
fprintf('Setup with time : %d s.\n', t_setup);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE - choose mmv or individual spg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don't compute validation error. (Used for PC calibration)
pc_val = 0; 

%%% reference solution
[c_ref, psi_ref, ~] = my_pce(xi_ref, p, u_ref, sigma, pc_solver, pc_val); 

t_ref = toc - t_setup; 
fprintf('Reference solution : %d s.\n', t_ref);


%%% Low-fidelity solution

[c_low, psi_low, ~] = my_pce(xi_low, p, u_low, sigma, pc_solver, pc_val); 

t_low = toc - t_ref; 
fprintf('Low fidelity solution : %d s.\n', t_low);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relative errors reference and low
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mean of reference solution
mean_ref = c_ref(:,1); 

% mean of low fidelity solution
mean_low = c_low(:,1);

var_ref = sum(c_ref(:,2:end).^2,2); 
std_ref = sqrt(sum(c_ref(:,2:end).^2,2)); 

% variance of low fidelity solution
var_low = sum(c_low(:,2:end).^2,2); 
std_low = sqrt(sum(c_low(:,2:end).^2,2)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bi-fidelity 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_points = size(u_ref,2); 	
n_samps = size(u_ref,1); 

P = size(c_ref,2); 

% Nodal covariance of the reference model
cov_ref = c_ref(:, 2:end)*c_ref(:,2:end)';
% Eigenvalue decomposition of nodal covariance
lam_ref = eigs(cov_ref, length(cov_ref));

sample = datasample(1:n_samps, N_hi, 'Replace', false); 

%%% High fidelity model - limited samples N_hi
xi_hi = xi_ref(sample,:); 
u_hi = u_ref(sample,:); 
psi_hi= psi_ref(sample,:); 

opts = spgSetParms('iterations',10000,'verbosity',0);

if pc_solver == 0
    c_hi = psi_hi\u_hi; 
    c_hi = c_hi';
elseif pc_solver == 1    
    % Solve PCE coefficents via l_1 minimization
    c_hi = spg_mmv(psi_hi,u_hi,sigma*norm(u_hi),opts);            
    c_hi = c_hi';
else
        weights = get_matrix_weights(psi_hi);
        Psiweights = psi_hi*weights;

        norm_u_vec = vecnorm(u_hi,2); 
%                     c_data{n_points} = []; 

        c_hi = zeros(n_points, P);
        for i_points=1:n_points
            opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

            delta = sigma*norm_u_vec(i_points);
            c_hi(i_points,:)  = weights*spg_bpdn(Psiweights,u_hi(:,i_points),delta,opts);
        end
end 

% size c_low is n_grid points (9) by basis functions P
cov_low = c_low(:, 2:end) * c_low(:, 2:end)'; 

% Projection of low-fid PC solution onto the eigenvectors of covariance
% yields size r (approximation rank) by P -1 
[v_low, lam_low] = eigs(cov_low, r); 

alpha2 = diag(1./sqrt(diag(lam_low)))*(v_low'*c_low(:, 2:end)); 

% New reduced basis (note: psi_hi is independent of low or high fid. Just
% depends on number of samples available. 
psi_bi = psi_hi(:,2:end)*alpha2';

% New reduced basis including column of ones
psi_bi = [ones(size(u_hi, 1), 1) psi_bi]; 


opts = spgSetParms('iterations', 10000, 'verbosity', 0); 

c_bi = spg_mmv(psi_bi, u_hi, sigma*norm(u_hi), opts); 
% c_bi = psi_bi\u_hi; 
c_bi = c_bi'; 

% Compute mean, variance or bi-fid estimate

% psi_ref and psi_low are the same. 

psi_bi_est = psi_ref(:,2:end)*alpha2';
psi_bi_est = [ones(size(u_ref, 1), 1) psi_bi_est]; 

u_bi = psi_bi_est*c_bi';

save('Results/LDC_qoi_mean_var', 'x_h', 'x_l', 'u_ref', 'u_low', 'u_bi', 'u_hi')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Change axes, but first update low-fidelity data... 
% % Then put everthing in slides to determine what to include. 
% % Then do low-fidelity bound check. 
% 
% 
% % biggest improvement
% [~, i_plot] = max(vecnorm(u_ref- u_low,2,2)-vecnorm(u_ref- u_bi,2,2));
% 
% 
% figure
% p1 = plot(x_h, u_ref(i_plot,:),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(x_l, u_low(i_plot,:),'-','color',c2,'LineWidth',LW);
% p3 = plot(x_h, u_bi(i_plot,:),'--','color',c3,'LineWidth',LW);
% xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
% ylabel('$u$','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% set(gcf, 'Position', size_1)
% legend([p1,p2,p3],{'Ref','L', 'B'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% if save_on == 1
%     saveas(gcf,'Plots/LDC_realization','png')
% end
% 
% % I need to check these coordinates... 
% 
% % Better to plot a realization, or the mean and the variance? 
% 
% % plot both for each and then decide. 
% 
% figure
% subplot(1,2,1)
% p1 = plot(x_h, mean(u_ref),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(x_h, mean(u_low),'-','color',c2,'LineWidth',LW);
% p3 = plot(x_h, mean(u_bi),'--','color',c3,'LineWidth',LW);
% xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
% ylabel('$u$ mean','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% 
% subplot(1,2,2)
% p1 = plot(x_h, var(u_ref),'-','color',c1,'LineWidth',LW);
% hold on
% p2 = plot(x_h, var(u_low),'-','color',c2,'LineWidth',LW);
% p3 = plot(x_h, var(u_bi),'--','color',c3,'LineWidth',LW);
% xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
% ylabel('$u$ variance','interpreter', 'latex', 'fontsize', FS)
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% legend([p1,p2,p3],{'Ref','L', 'B'},'interpreter', 'latex', 'fontsize', FS_leg)
% 
% set(gcf, 'Position', size_2)
% 
% if save_on == 1
%     saveas(gcf,'Plots/LDC_mean_var','png')
% end


