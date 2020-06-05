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

% Presentation size
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

save_on = 1; 

QoI = 0; % u mid
% QoI = 1; % cylinder


if QoI == 0
    results_name = 'GT_mid_';
    label_name = 'Vertical Line';
elseif QoI == 1
    results_name = 'GT_cylinder_'; 
    label_name = 'Cylinder Surface';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load data                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if QoI == 0
    load('Results/GT_u_mid_results_spg.mat')
elseif QoI == 1
    load('Results/GT_cylinder_results_spg.mat')
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

ylabel('Average Relative Error in Mean','interpreter','latex','Fontsize',FS)
xlabel('$n$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
set(gcf, 'Position', size_1)
% grid on
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
ylabel('Average Relative Error in Variance','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)
% grid on 
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
xlabel('Index $i$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Normalized Eigenvalue', 'interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([1,10])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2, p3, p4],{'L QoI 1','Ref QoI 1', 'L QoI 2','Ref QoI 2'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'eigen'),'epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QoI mean and variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat('Results/',results_name, 'mean_var'));

figure
subplot(1,2,1)
p0 = plot(x_h, mean(u_ref),'k:+','LineWidth',LW);
hold on
p1 = plot(x_h, mean(u_hi),'-','color',c1,'LineWidth',LW);
p2 = plot(x_l, mean(u_low),'--','color',c2,'LineWidth',LW);
p3 = plot(x_h, mean(u_bi),'-.','color',c3,'LineWidth',LW);
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
p0 = plot(x_h, var(u_ref),'k:+','LineWidth',LW);
hold on
p1 = plot(x_h, var(u_hi),'-','color',c1,'LineWidth',LW);
p2 = plot(x_l, var(u_low),'--','color',c2,'LineWidth',LW);
p3 = plot(x_h, var(u_bi),'-.','color',c3,'LineWidth',LW);
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
legend([p0,p1,p2,p3],{'Ref','H','L', 'B'},'interpreter', 'latex', 'fontsize', FS_leg)
% title('Variance','interpreter', 'latex', 'fontsize', FS_axis)

set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'mean_var'),'epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound efficacy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/GT_efficacy_cyl');
efficacy_cy = efficacy; 
% efficacy_cy = efficacy_mid; % Until fixed

load('Results/GT_efficacy_mid');
efficacy_mid = efficacy; 

% Set common limits 
lim_min = min([efficacy_cy(:); efficacy_mid(:)]);
lim_max = max([efficacy_cy(:); efficacy_mid(:)]);


figure
subplot(1,2,1)
h = pcolor(N_hi_vec, r_vec, efficacy_mid);
set(h, 'EdgeColor', 'none');
axis tight
xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Approximation rank $r$', 'interpreter', 'latex', 'fontsize', FS)
colorbar
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
title('Error Bound Efficacy: Vertical Line','interpreter', 'latex', 'fontsize', FS_axis)
caxis([lim_min, lim_max]); 

subplot(1,2,2)
h = pcolor(N_hi_vec, r_vec, efficacy_cy);
set(h, 'EdgeColor', 'none');
axis tight
xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Approximation rank $r$', 'interpreter', 'latex', 'fontsize', FS)
colorbar
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
caxis([lim_min, lim_max]); 
title('Error Bound Efficacy: Cylinder Surface','interpreter', 'latex', 'fontsize', FS_axis)

set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,'Plots/GT_efficacy','epsc')
end
