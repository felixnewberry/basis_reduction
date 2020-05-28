%%% Test sampling bound

close all
clear all
clc

%%%

% Estimating airfoil data. 
% 128 grid points, N = 1200

% Computed A_n for n = 1000. 
% Psi_bi found with p = 5, r = 30. Should have minimal sampling error. 


% At present q < 0. 
% Intuition suggests this should be a sufficient number of samples. 

% Possible reason is that q computes a bound for an individual estimate,
% not a matrix estimate. 

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
%%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 1200; 

load('A.mat')

% Bi-fidelity basis found with p = 5, r = 30. 
% Should have minimal sampling error. 
% load('psi_bi_n1000_r30_N.mat')
% psi_bi_1000 = psi_bi_N; 
% 
% load('psi_bi_n200_r30_N.mat')
% psi_bi_200 = psi_bi_N; 

load('psi_bi_n1000_r30_n.mat')
psi_bi_1000 = psi_bi_n; 

load('psi_bi_n200_r30_n.mat')
psi_bi_200 = psi_bi_n; 

% First n hi-fidelity samples. 128 grid points
load('A_n1000.mat')
A_n1000 = A_n; 

load('A_n200.mat')
A_n200 = A_n; 

% load bi-fidelity estimate
load('bi_est_n1000_r30')

% n = 200 or 1000

% %%% Should be close to 0
% psi_bi = psi_bi_1000; 
% A_n = A_n1000; 
% n = 1000; 

%%% Should be well above 0
psi_bi = psi_bi_200; 
A_n = A_n200; 
n = 200; 

% If I change mu from being normalized to not it goes from impossibly tight
% to far too loose. 

% compute psi_bi anew and see what happens. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute point bound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% omega is an upper bound on realizations of u. 
% compute with surrogate of subset of n high-fidelity samples. 

% omega - maximum of each sample - estimate with H_n as surrogate
omega_term = sum(max(A_n,[],2).^2)^0.5;

% A1 = mat_normalize(psi_bi(:,2:end)');
A1 = psi_bi(:,2:end); 



mu = max(sum(abs(A1).^2,2)); 

%
M = 1/size(A1,1)*(A1'*A1); 
s = norm(M - eye(length(M)));

% > 1/2 - does not satisfy condition coherence motivated sampling. 
%
% Lemma 1: 
size(A1,1)*sqrt(size(A1,2)); %>= 20
Prob = 2*size(A1,2)*exp(-0.1*size(A1,1)*mu^-1);


% check univariate vs multivariate basis in paper with Jarred
% is the ratio from cohen's paper still valid? 

mu_normal = max(sum(abs(mat_normalize(psi_bi(:,2:end)')).^2,1)); 

q_N = 0.5*(N*(3*log(3/2)-1)/(mu*log(N))-2); 
q_n = 0.5*(n*(3*log(3/2)-1)/(mu*log(n))-2);

% I should put a check here on sign
if q_N <= 0
    fprintf("Error: For N = %d, q_N < 0, (q_n = %6.2f), and condition is invalid. \n", N, q_N);
end
if q_n <= 0
    fprintf("Error: For n = %d, q_N < 0, (q_n = %6.2f), and condition is invalid. \n", n, q_n);
end

err_bound_N = 2*sqrt(2)*N^0.5*omega_term*(N^(-q_N/2))
err_bound_n = 2*sqrt(2)*N^0.5*omega_term*(n^(-q_n/2))


%%% bound of just one point: 

omega_points = abs(max(A_n,[],2)); 

bound_points = 2*sqrt(2)*omega_points*n^(-q_n/2); 
% Plot against the points themselves

figure
p1 = plot(A,'-','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(A_hat,'--','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = plot(bound_points,'.-','Color',c3,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('points', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
legend([p1(1),p2(1), p3(1)],{'H','B', 'Samp Bound',},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

% Looks very loose... 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Vary n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If I normalize omega the bound is too small, but q is positive. 
% Mu = 1 when bound is normal. 

% If I don't normalize omega, the bound is massive - and gets bigger as N
% increases because q <0 

n_vec = 50:50:1000; 
err_n_vec = zeros(1,length(n_vec)); 
q_n_vec = zeros(1,length(n_vec)); 

for i_n = 1:length(n_vec)
    n = n_vec(i_n); 
    
    A_n = A(:,1:n); 
    
    % omega - maximum of each sample - estimate with H_n as surrogate
    omega_term = sum(max(A_n,[],2).^2)^0.5;
    
%     A1 = mat_normalize(psi_bi(:,2:end)');
    A1 = psi_bi(:,2:end)'; 
    mu = max(sum(abs(A1).^2,1)); 

    q_n_vec(i_n) = 0.5*(n*(3*log(3/2)-1)/(mu*log(n))-2);
    err_n_vec(i_n) = 2*sqrt(2)*N^0.5*omega_term*(n^(-q_n_vec(i_n)/2));
end

figure
p1 = semilogy(n_vec, err_n_vec,'-','Color',c1,'LineWidth',LW,'MarkerSize',MS);
xlabel('$n$ samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error Bound', 'interpreter', 'latex', 'fontsize', FS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A LS problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
