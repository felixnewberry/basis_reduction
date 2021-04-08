clear all 
close all
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot Settings                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
MS = 8;     % Marker Size
FS_leg = 16; % Font size legend

% % Presentation size
% size_1 = [0,0,445,345]; 
% size_2 = [0,0,890,345]; 

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

save_on = 0; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Choose grid points to look at, 128 total.
% gridpt = [1:128];
% 
% load('x_locations')
% x_l = x_locations; 
% x_h = x_locations;
% 
% load lowFiResults_6d2.mat;
% u_low = lowFiResults(gridpt,:)'; % Number of samples x spatial/temporal dimension of problem
% 
% load highFiResults_6d2.mat;
% u_ref = highFiResults(gridpt,:)';
% 
% % load random variables used to generate u
% load Inputs_6d2.mat;
% 
% %     normalize to be U~[-1,1] and re-label as xi.
% % xi_ref should be of dimension: Number of reference samples x stochastic
% % dimension d
% xi_ref =[(1/3)*alpha',(4*c-4)',(50*m-1)',(10*(mach-.2)-2)',(10*p -4)',((t - .1)*50 - 1)'];
% xi_low = xi_ref; % number of low fidelity samples X stochastic dimension d
% 
% 
% %%% Interpolate data to be evenly distributed
% % Adjust coordinates to go from -1 to 1. 
% % clockwise from trailing edge: positive is pressure surface (lower)
% 
% x_order = [1-x_h(65:end);x_h(1:64)-1];
% u_low_order = [u_low(:,65:end), u_low(:,1:64)];
% u_ref_order = [u_ref(:,65:end), u_ref(:,1:64)];
% 
% % u_low and u_ref
% n_points = 100; 
% n_sim = length(u_ref); 
% 
% x_int = linspace(-1,1,n_points); 
% u_low = zeros(n_sim, n_points); 
% u_ref = zeros(n_sim, n_points); 
% for i_sim = 1:n_sim
%     u_low(i_sim,:) = interp1(x_order, u_low_order(i_sim,:), x_int);
%     u_ref(i_sim,:) = interp1(x_order, u_ref_order(i_sim,:), x_int);
% end

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

N_hi = 20;      % Number high-fidelity samples
% N_hi = [30]; 

r = 8;           % KL order

% tolerance on residual used in spgl1 solver
sigma = .001;

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

% mean_low_int = interp1q(x_l', mean_low, x_h');
mean_low_err = norm(mean_low - mean_ref)/norm(mean_ref); 

% var_low_int = interp1q(x_l', var_low, x_h');
var_low_err = norm(var_low - var_ref)/norm(var_ref); 

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
psi_bi = [ones(size(u_hi, 1), 1) psi_bi(:,1:end-1)]; 


opts = spgSetParms('iterations', 10000, 'verbosity', 0); 

c_bi = spg_mmv(psi_bi, u_hi, sigma*norm(u_hi), opts); 
% c_bi = psi_bi\u_hi; 
c_bi = c_bi'; 

% Compute mean, variance or bi-fid estimate


save('Results/Airfoil_qoi_mean_var', 'x_int', 'c_ref', 'c_low', 'c_bi', 'c_hi')


figure
subplot(1,2,1)
p0 = plot(x_int, c_ref(:,1),'k:+','LineWidth',LW);
hold on
p1 = plot(x_int, c_hi(:,1),'-','color',c1,'LineWidth',LW);
p2 = plot(x_int, c_low(:,1),'--','color',c2,'LineWidth',LW);
p3 = plot(x_int, c_bi(:,1),'-.','color',c3,'LineWidth',LW);
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$ Mean','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on

subplot(1,2,2)
p0 = plot(x_int, sum(c_ref(:,2:end).^2,2),'k:+','LineWidth',LW);
hold on
p1 = plot(x_int, sum(c_hi(:,2:end).^2,2),'-','color',c1,'LineWidth',LW);
p2 = plot(x_int, sum(c_low(:,2:end).^2,2),'--','color',c2,'LineWidth',LW);
p3 = plot(x_int, sum(c_bi(:,2:end).^2,2),'-.','color',c3,'LineWidth',LW);
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$ Variance','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p0,p1,p2,p3],{'Ref','$H$', '$L$', '$B$'}, 'interpreter', 'latex', 'fontsize', FS_leg)

set(gcf, 'Position', size_2)