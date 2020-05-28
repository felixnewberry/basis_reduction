clear all 
close all 
clc

% Plot Gas Turbine

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

save_on = 0; 

QoI = 1; % u mid
% QoI = 1; % cylinder


if QoI == 0
    results_name = 'GT_u_mid_';
elseif QoI == 1
    results_name = 'GT_cylinder_'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load data                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if QoI == 0
    load('Results/GT_u_mid_results.mat')
elseif QoI == 1
    load('Results/GT_cylinder_results.mat')
end

% Vector of strings for r plots
r_symbol = {'-.+','-.*','-.s','-.d'}; 

r_string = cell(length(r)+2,1);
r_string(1:2,:) = {'$H$'; '$L$'};
r_string(3:end) = cellstr(num2str(r', 'B(r=%-d)')); 

r_string_variance = r_string([1,3:end]); 

% r_plot = [1 2 3];
i_eigen = length(r);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mean and Variance            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bi_stats = cat(3,mean_mean_bi_err , mean_mean_hi_err, var_mean_bi_err, ...
%     var_mean_hi_err, mean_var_bi_err , mean_var_hi_err ,...
%     var_var_bi_err, var_var_hi_err);
% Plot together

mean_mean_hi_err = bi_stats(:,:,2);
mean_mean_bi_err = bi_stats(:,:,1);

mean_var_hi_err = bi_stats(:,:,6);
mean_var_bi_err = bi_stats(:,:,5);


figure

subplot(1,2,1)
semilogy(N_hi,(mean_mean_hi_err(1,:)'),'-o','Color',c1, ...
    'LineWidth',LW,'MarkerSize',MS)
hold on
semilogy(N_hi,repmat(mean_low_err, 1, length(N_hi)),'--x','Color',c2, ...
    'LineWidth',LW,'MarkerSize',MS)
hold on
% semilogy(N_hi,(mean_mean_bi_err'),'-ob','LineWidth',LW)
for i_r = 1:length(r)
    semilogy(N_hi,(mean_mean_bi_err(i_r,:)),...
        r_symbol{i_r},'Color',c3, 'LineWidth',LW,'MarkerSize',MS)
end

ylabel('Average relative errors in mean','Fontsize',FS)
xlabel('$n$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_1)
grid on
set(gcf, 'Position', size_1)

subplot(1,2,2)
semilogy(N_hi,(mean_var_hi_err(1,:)'),'-o','Color',c1, ...
    'LineWidth',LW,'MarkerSize',MS)
hold on
semilogy(N_hi,repmat(var_low_err,1,length(N_hi)),'--*','Color',c2, 'LineWidth',LW)
hold on

for i_r = 1:length(r)
    semilogy(N_hi,(mean_var_bi_err(i_r,:)),r_symbol{i_r},...
        'Color',c3, 'LineWidth',LW,'MarkerSize',MS)
end
xlabel('$n$','interpreter','latex','Fontsize',FS)
ylabel('Average relative errors in variance','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)
grid on 
set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'N_r'),'epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Eigenvalues             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot both QoI on one figure

% QoI 0
load('Results/GT_u_mid_results.mat')
mean_lam_low_0 = mean_lam_low; 
mean_lam_ref_0 = mean_lam_ref; 

% QoI 1
load('Results/GT_cylinder_results.mat')
mean_lam_low_1 = mean_lam_low; 
mean_lam_ref_1 = mean_lam_ref; 

% eigenvalues recorded for for each r N_hi combination.
figure
p1 = semilogy(abs(mean_lam_low_0)./max(mean_lam_low_0),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(abs(mean_lam_ref_0)./max(mean_lam_ref_0),'--x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(abs(mean_lam_low_1)./max(mean_lam_low_1),'-s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p4 = semilogy(abs(mean_lam_ref_1)./max(mean_lam_ref_1),'-->','Color',c4,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\lambda_i$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
grid on
set(gcf,'Position',size_1)

legend([p1,p2, p3, p4],{'L QoI 1','Ref QoI 1', 'L QoI 2','Ref QoI 2'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'eigen'),'epsc')
end