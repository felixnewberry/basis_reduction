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
%%% Chose QoI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QoI = 0; % u mid
QoI = 1; % cylinder

if QoI == 0
    results_name = 'GT_mid_';
    label_name = 'Vertical Line';
elseif QoI == 1
    results_name = 'GT_cylinder_'; 
    label_name = 'Cylinder Surface';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if QoI == 0
    results_name = 'GT_mid_'; 
    lowFiResults = importdata('assembledRunMid_40');   
    lowFiResults = lowFiResults';

    highFiResults = importdata('assembledRunMid_110_high'); 
    highFiResults = highFiResults';

    x_h = importdata('assembledCoordsMid_110_high');
    [~,idx_h] = sort(x_h(2,:));
    x_h = x_h(2,idx_h); 
    highFiResults = highFiResults(idx_h,:); 
  
    x_l = importdata('assembledCoordsMid_40');
    [~,idx_l] = sort(x_l(2,:));
    x_l = x_l(2,idx_l);
    lowFiResults = lowFiResults(idx_l,:); 
    
elseif QoI == 1
    results_name = 'GT_cylinder_'; 
    lowFiResults = importdata('assembledRunCylinder_40'); 
    lowFiResults = lowFiResults';
    tic

    highFiResults = importdata('assembledRunCylinder_110_high'); 
    highFiResults = highFiResults';
    
    x_h = importdata('assembledCoordsCylinder_110_high');
    x_h = atan2(x_h(2,:),-x_h(1,:));
    [~,idx_h] = unique(x_h);
    x_h = x_h(idx_h); 
    highFiResults = highFiResults(idx_h,:); 
  
    x_l = importdata('assembledCoordsCylinder_40');
    x_l = atan2(x_l(2,:),-x_l(1,:));    
    [~,idx_l] = unique(x_l);
    x_l = x_l(idx_l);
    lowFiResults = lowFiResults(idx_l,:);  
end

% y_samp = load('uniform_40_low.mat');
% xi_low = y_samp.uniform_data;

y_samp = load('uniform_110_high.mat');
xi_ref = y_samp.uniform_data;
xi_low = xi_ref; 
   

gridpt_l = 1:length(x_l); 
gridpt_h = 1:length(x_h); 

% gridpt_h = 1:4; 

u_ref = highFiResults(gridpt_h,:)';
u_low = lowFiResults(gridpt_l,:)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Key Parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 6;                          % PCE order
% d = 4;                          % Stochastic dimension

N_hi = 50;      % Number high-fidelity samples
% N_hi = [30]; 

r = 8;           % KL order

% tolerance on residual used in spgl1 solver
sigma = .007;

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
% mean_low_err = norm(mean_low - mean_ref)/norm(mean_ref); 

% var_low_int = interp1q(x_l', var_low, x_h');
% var_low_err = norm(var_low - var_ref)/norm(var_ref); 

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

% psi_ref and psi_low are the same. 

psi_bi_est = psi_ref(:,2:end)*alpha2';
psi_bi_est = [ones(size(u_ref, 1), 1) psi_bi_est(:,1:end-1)]; 

u_bi = psi_bi_est*c_bi';

% Check errors: 


mean_low = c_low(:,1);
mean_low_int = interp1q(x_l', mean_low, x_h');
mean_low_err = norm(mean_low_int - mean_ref)/norm(mean_ref); 

save(strcat('Results/', results_name, 'mean_var'), 'x_h', 'x_l', 'c_ref', 'c_low', 'c_bi', 'c_hi')

figure
subplot(1,2,1)
p0 = plot(x_h, c_ref(:,1),'k:+','LineWidth',LW);
hold on
p1 = plot(x_h, c_hi(:,1),'-','color',c1,'LineWidth',LW);
p2 = plot(x_l, c_low(:,1),'--','color',c2,'LineWidth',LW);
p3 = plot(x_h, c_bi(:,1),'-.','color',c3,'LineWidth',LW);
if QoI == 0
    xlabel('$y$','interpreter', 'latex', 'fontsize', FS)
else
    xlabel('$\theta$','interpreter', 'latex', 'fontsize', FS)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    xlim([-pi, pi])
end
ylabel(strcat(label_name, ' Temperature Mean'),'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% title('Mean','interpreter', 'latex', 'fontsize', FS_axis)

subplot(1,2,2)
p0 = plot(x_h, sum(c_ref(:,2:end).^2,2),'k:+','LineWidth',LW);
hold on
p1 = plot(x_h, sum(c_hi(:,2:end).^2,2),'-','color',c1,'LineWidth',LW);
p2 = plot(x_l, sum(c_low(:,2:end).^2,2),'--','color',c2,'LineWidth',LW);
p3 = plot(x_h, sum(c_bi(:,2:end).^2,2),'-.','color',c3,'LineWidth',LW);
if QoI == 0
    xlabel('$y$','interpreter', 'latex', 'fontsize', FS)
else
    xlabel('$\theta$','interpreter', 'latex', 'fontsize', FS)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'})
    xlim([-pi, pi])
end
ylabel(strcat(label_name, ' Temperature Variance'),'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p0,p1,p2,p3],{'Ref','$H$','$L$', '$B$'},'interpreter', 'latex', 'fontsize', FS_leg)
% title('Variance','interpreter', 'latex', 'fontsize', FS_axis)

set(gcf, 'Position', size_2)
