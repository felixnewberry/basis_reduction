clear all 
close all 
clc

% Plot Airfoil

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

save_on = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load data                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Results/Airfoil_results.mat')
% load('Results/Airfoil_results_spg.mat')
% load('Results/Airfoil_results_spg_int.mat')
load('Results/Airfoil_results_spg_int_2.mat')

% load('Results/Airfoil_results_int.mat')

% Vector of strings for r plots
r_symbol = {'-.+','-.*','-.s','-.d'}; 

r_string = cell(length(r)+2,1);
r_string(1:2,:) = {'$H$'; '$L$'};
r_string(3:end) = cellstr(num2str(r', 'B(r=%-d)')); 

r_string_variance = r_string([1,3:end]); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Eigenvalues             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% abs to account for machine error close to 0.
figure
p1 = semilogy(abs(mean_lam_low)./max(mean_lam_low),'-o','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(abs(mean_lam_ref)./max(mean_lam_ref),'--x','Color',c2,...
    'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('Index $i$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Normalized Eigenvalue', 'interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([1,10])
yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2],{'L','Ref'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,'Plots/Airfoil_eigen','epsc')
end

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
    saveas(gcf,'Plots/Airfoil_N_r','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Airfoil QoI Mean and Variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/Airfoil_qoi_mean_var')

figure
subplot(1,2,1)
p0 = plot(x_int, mean(u_ref),'k:+','LineWidth',LW);
hold on
p1 = plot(x_int, mean(u_hi),'-','color',c1,'LineWidth',LW);
p2 = plot(x_int, mean(u_low),'--','color',c2,'LineWidth',LW);
p3 = plot(x_int, mean(u_bi),'-.','color',c3,'LineWidth',LW);
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$ mean','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on

subplot(1,2,2)
p0 = plot(x_int, var(u_ref),'k:+','LineWidth',LW);
hold on
p1 = plot(x_int, var(u_hi),'-','color',c1,'LineWidth',LW);
p2 = plot(x_int, var(u_low),'--','color',c2,'LineWidth',LW);
p3 = plot(x_int, var(u_bi),'-.','color',c3,'LineWidth',LW);
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$ variance','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend([p0,p1,p2,p3],{'Ref','H', 'L', 'B'}, 'interpreter', 'latex', 'fontsize', FS_leg)

set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,'Plots/Airfoil_mean_var','epsc')
end

% I could reverse the axis with  'YDir', 'reverse' in set(gca.. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Airfoil Efficacy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


load('Results/Airfoil_efficacy');

figure
h = pcolor(N_hi_vec, r_vec, efficacy);
set(h, 'EdgeColor', 'none');
axis tight
xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Approximation rank $r$', 'interpreter', 'latex', 'fontsize', FS)
colorbar
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% title('Efficacy','interpreter', 'latex', 'fontsize', FS_leg)
set(gcf, 'Position', size_1)

if save_on == 1
    saveas(gcf,'Plots/Airfoil_efficacy','epsc')
end

