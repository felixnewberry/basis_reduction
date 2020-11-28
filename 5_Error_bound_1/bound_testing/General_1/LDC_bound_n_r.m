% Test Error Bound

% Problem setup

clear all
close all
clc

tic 

save_on = 1; 

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

results_file_name = 'LDC_bound'; 

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

%%% reference solution
[c_ref, psi_ref] = my_pce(xi_ref, p, A', sigma, pc_solver); 
t_ref = toc(timerval); 
fprintf('Reference solution : %d s.\n', t_ref);
timerval = tic; 

%%% Low-fidelity solution
[c_low, psi_low] = my_pce(xi_low, p, B', sigma, pc_solver); 
t_low = toc(timerval); 
% t1 = toc
% t_ref
fprintf('Low fidelity solution : %d s.\n', t_low);


%%% Low_fidelity 
u_low_int = interp1(x_l, B, x_h);  
err_low = norm(u_low_int-A)/norm(A); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over N_hi, or n to compute epsilon tau, or r for reduced basis... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I want a function that does these things. 

% N_hi_vec = 35:5:40; 
% N_hi = 80; 
N_hi = 40; 

n_vec = N_hi-10:N_hi+10;

%%% Check that q > 0 
% works for 33 for LDC - Bound is very loose at this point. 
% Doesn't seem useful in general... 

% N_hi_vec = 3:4; 
% r = 3; 
r_vec = 3:7; 

length_n = length(n_vec);
length_r = length(r_vec); 

% n_reps = 1; 

err_bi_vec = zeros(length_r,length_n); 
err_bound_vec = zeros(length_r,length_n); 
err_ahat_vec = zeros(length_r,length_n);  
err_bound_N_vec = zeros(length_r,length_n); 
err_bound_n_vec = zeros(length_r,length_n); 

for i_n = 1:length_n
    for i_r = 1:length_r

    [err_bi_vec(i_r, i_n), err_bound_vec(i_r, i_n), err_ahat_vec(i_r, i_n), ...
        err_bound_N_vec(i_r, i_n), err_bound_n_vec(i_r, i_n), ~, ...
        ~, ~, ~] = ...
        br_bi_bound(B, A, N_hi, n_vec(i_n), psi_ref, psi_low, c_low, sigma, ...
        r_vec(i_r), p, xi_low, pc_solver); 
    
    end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Efficacy 

efficacy = err_bound_vec./err_bi_vec

% 3 when 80 N_hi samples used. Little variation. 
% 84 when 40 used: Predicting 7 %... 

% Try other problems and see if behaviour is consistent. 

% High efficacy - but requires loads of samples :( 
% %%% Bi and bound
% figure
% p1 = semilogy(N_hi_vec,err_bi_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
% hold on
% p2 = semilogy(N_hi_vec,err_bound_vec,'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
% hold off
% axis tight
% xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% % grid on
% set(gcf,'Position',size_1)
% legend([p1,p2],{'Bi','Bound'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
% title('Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)
% 
% if save_on == 1
%     saveas(gcf,'bound_plots/LDC_bound_performance','png')
% end    
