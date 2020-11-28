%%% univariate test case

clear all
close all
clc

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
%%% Univariate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(42); 

N_total = 1e5; 

N_train = N_total; 

d = 1; 
p_ref = 3; 

xi_ref = rand(N_train,d)*2-1;

%%% Set up data
% c_ref = [1, 0.4, 0.1, 0.2]'; 
c_ref = [1, 0.4, 0.1, 0.02]'; 



% index_pc = nD_polynomial_array(d,p); 
index_pc_ref = [0:p_ref]';

% size of PC basis (set of P basis functions 
P_ref = size(index_pc_ref,1);

% Construct reference polynomial basis
psi_ref = zeros(N_train,P_ref);    
for isim=1:N_train
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_ref(isim,:),index_pc_ref);
    psi_ref(isim,:) = crow_ref(1:P_ref);
end

u_ref = psi_ref*c_ref; 


%%% estimate with p = 2... 

p_est = 2; 
% index_pc = nD_polynomial_array(d,p); 
index_pc_est = [0:p_est]';

% size of PC basis (set of P basis functions 
P_est = size(index_pc_est,1);

% Construct reference polynomial basis
psi_est = zeros(N_train,P_est);   

for isim=1:N_train
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_est = piset(xi_ref(isim,:),index_pc_est);
    psi_est(isim,:) = crow_est(1:P_est);
end

% Decide if using correct or incorrect model: 
% psi_use = psi_ref; 
psi_use = psi_est; 

% Error is massive... I want to simulate truncation error. make the higher
% order term smaller? 

% Still - bound should work - fix up code to do so. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loop over error and bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_reps = 100; 
sigma = 1e-3;
% n_samp_vec = 15:5:50; 
n_samp_vec = 4:1:7;  % r_opt = 5.2
n_samp_vec = 4:1:30;  % r_opt = 5.2

opts = spgSetParms('iterations',10000,'verbosity',0);

% store squared errors
err_n_mat = zeros(n_reps,length(n_samp_vec)); 

% Store terms 1 and 2 of cohen et al. bound
% bound_cohen_t1 = zeros(n_reps,length(n_samp_vec)); 
bound_cohen_t2 = zeros(n_reps,length(n_samp_vec)); 

% Store coherence motivated sampling bound
bound_pc_mat =  zeros(n_reps,length(n_samp_vec)); 

for i_rep = 1:n_reps
    for i_n = 1:length(n_samp_vec)
        
        %%% Approximation
        % Select sample
        n_samps = n_samp_vec(i_n); 
        rand_sample = datasample(1:N_train,n_samps,'Replace',false); 
        u = u_ref(rand_sample,1);
        xi = xi_ref(rand_sample); 
        

        % solve for coefficients with sample via least squares
        psi = psi_use(rand_sample,:); 
        c = psi\u; 

        % Estimate the training data
        u_est = psi_use*c; 

        err_n_mat(i_rep, i_n) = mean((u_est - u_ref).^2); 
        
        % Can use coefficients if the basis is correct
%         err_n_mat(i_rep, i_n) = norm(c-c_ref).^2; 

        %%% Cohen et al bound: 
        
        % Coherence:
        A1 = psi(:,2:end);  % should I be doing this? 
        
        mu_cohen = max(sum(abs(A1).^2,2)); 
        mu_normal = max(sum(abs(mat_normalize(psi(:,2:end)')).^2,1)); 
        mu_cohen = mu_normal; 
        
        % Bound not valid if q_n <0. Happens if mu is not normalized
        q_n = 0.5*(n_samps*(3*log(3/2)-1)/(mu_cohen*log(n_samps))-2);
        
        omega = max(abs(u)); 
        bound_cohen_t2(i_rep, i_n)=  8*omega^2*n_samps^(-q_n); 
        
        % best approximation error is 0
 
        % q is still negative (if coherence not normalized). What am I doing wrong?? 
        
        %%% Coherence motivated sampling bound
        A1 = psi(:,2:end);  % should I be doing this?
        mu_pc = max(sum(abs(A1).^2,2)); 
        
        % Surrogate truncation error? Think about this... 
        bound_pc_mat(i_rep, i_n) = err_n_mat(i_rep, i_n) *(1 + 4*mu_pc./n_samps); 
    end

end

%%% Check M condition: - this wasn't satisfied in the multivariate case
M = 1/size(psi_ref,1)*psi(:,2:end)'*psi(:,2:end); 
s = norm(M - eye(size(M))); 

% s < 0.5 - I can proceed with bounds. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% calculate mean and variances
err_n_vec = mean(err_n_mat,1); 
err_n_var = var(err_n_mat,1); 

bound_vec_t2 = mean(bound_cohen_t2,1); 
bound_pc = mean(bound_pc_mat,1); 


% n_use for r
n_r_calc = 5; 

% Find r through gradient: 
coefs = polyfit(n_samp_vec(1:n_r_calc), log(err_n_vec(1:n_r_calc)), 1);

figure
p1 = plot(n_samp_vec(1:n_r_calc), log(err_n_vec(1:n_r_calc)),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = plot(n_samp_vec(1:n_r_calc), coefs(2)+coefs(1)*n_samp_vec(1:n_r_calc),'--','Color',c2,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\log{(\Vert e \Vert^2)}$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Optimal r calculation','interpreter', 'latex', 'fontsize', FS_leg)

% r_opt = -coefs(1)
% r_opt = 0.0672; 
% r_opt = 5.2366; % p = 3
% r_opt = 5.35; % p = 4

% Guess at optimal... 
% r_opt = 3; 
r_opt = 0.4; 

omega = max(abs(u_ref)); 
bound_cohen_t2=  8*omega^2*n_samp_vec.^(-r_opt);  

kappa = (3*log(3/2)-1)/(2+2*r_opt); 
epsilon = 4*kappa./log(n_samp_vec); 

e_mf2 = err_n_vec(end).^2; 
bound_cohen_t1 = (1+epsilon).*e_mf2; 

% Plot squared error
figure
p1 = semilogy(n_samp_vec, err_n_vec,'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(n_samp_vec, err_n_vec+sqrt(err_n_var),'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
semilogy(n_samp_vec, err_n_vec-sqrt(err_n_var),'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(n_samp_vec, bound_pc,'-.s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p4 = semilogy(n_samp_vec, bound_cohen_t1,'-.<','Color',c4,'LineWidth',LW,'MarkerSize',MS);
p5 = semilogy(n_samp_vec, bound_cohen_t2,'-.*','Color',c5,'LineWidth',LW,'MarkerSize',MS);
p6 = semilogy(n_samp_vec, bound_cohen_t1+bound_cohen_t2,'-.d','Color',c6,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('n samples', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error $^2$', 'interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p1,p2, p3, p4, p5, p6],{'Bi','$\sigma$', 'PC Bound', 'Cohen t1',  'Cohen t2', 'Cohen total'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
title('Sampling Bound Performance','interpreter', 'latex', 'fontsize', FS_leg)


%%% Try a more comprehensive fit... 

r_opt = 3; 
log_b = log(bound_cohen_t1+8*omega^2*n_samp_vec.^(-r_opt)); 
log_e = log(err_n_vec); 

figure
plot(log_e, '-bx')
hold on
plot(log_b, '-.ro')