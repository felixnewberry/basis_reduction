

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;     % Line width
FS = 24;    % Font size axis
FS_leg = 16; % Font size legend
FS_axis = 14; 
LW_axis = 2; 

FS = 28;    % Font size axis
FS_axis = 18; 
LW_axis = 2; 
size_1 = [0, 0, 410, 300]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem Specific
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ie random variables are UB_vel, LB_vel and Visc
% or UB_mean, UB_amplitude and UB_freq

% Dimension (d) and total order (p) of PC
d = 3; % 3 for UB_LB_visc
p = 3; % 3 for UB_LB_visc 

% Number of samples N

% UB_LB_visc
% N_total = 100;
% N = 70;
% nsim_v = N_total - N;


% UB_mean UB_amp UB_freq
N_total = 50;
N = 40;
nsim_v = N_total - N;


% Load the PCE idex matrix
index_pc = nD_polynomial_array(d,p);

% size of the basis
P = size(index_pc,1);

% Load data

load assembledData_UB_LB_visc;
u_samp = assembledData_UB_LB_visc(1:N_total,1);
load uniform_data_UB_LB_visc.mat;
y_samp = uniform_data(1:N_total,:); % samples of the input variables

% Solve the l1 minimization problem min || c ||_1 s.t. ||Psi c - u_samp||_2^2 < delta 

% Form the measurement matrix Ps
for isim=1:N
    crow = piset(y_samp(isim,:),index_pc);
    Psi(isim,:) = crow;
end

% Solve for the solution using least squares
c_ls = Psi\u_samp(1:N,1);

% Plot the coefficients 
loglog(abs(c_ls),'ro'); 
legend('l2 minimization')

% Relative mean square error
%err_rms = norm(c_l1-c_ref(1:P))/norm(c_ref(1:P)); disp(err_rms);


%% Start validation

clear Psi;
for isim=1:nsim_v
    crow = piset(y_samp(isim+N,:),index_pc);
    Psi(isim,:) = crow(1:P);
end
error_val = norm(u_samp(N+1:end,1)-Psi*c_ls)/norm(u_samp(N+1:end,1))

%% Comput the Sobol' indices
% Call the Sobol index code
for j=1:d    
    s(j) = sobol_pce( c_ls, index_pc(:,j));
end

% Plot Sobol' indices
labels = {' \xi_1 ','    \xi_2 ',' \xi_3    '};
explode = heaviside(s);
figure;
p = pie(s,explode,labels);
colormap([255 244 102; %// yellow
        112 193 129;      %// green
        36 123 160]./255) %// blue
t1 = p(2); 
t2 = p(4); 
t3 = p(6); 

t1.FontSize = 18;
t2.FontSize = 18; 
t3.FontSize = 18; 

% set(gca,'Fontsize', FS_axis);
% xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$\lambda_i$', 'interpreter', 'latex', 'fontsize', FS)
% set(gca,'linewidth',LW_axis);
% box on
set(gcf, 'Position', size_1)
% grid on
