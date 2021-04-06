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

save_on = 0; 

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
r_string(3:end) = cellstr(num2str(r', '$B(r=%-d)$')); 

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2],{'$L$','Ref'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)


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
yl = ylim; 
ylim([yl(1),1])
new_labels = [1e-1, 1];
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
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
p0 = plot(x_int, c_ref(:,1),'k:+','LineWidth',LW);
hold on
p1 = plot(x_int, c_hi(:,1),'-','color',c1,'LineWidth',LW);
p2 = plot(x_int, c_low(:,1),'--','color',c2,'LineWidth',LW);
p3 = plot(x_int, c_bi(:,1),'-.','color',c3,'LineWidth',LW);
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$ Mean','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([p0,p1,p2,p3],{'Ref','$H$', '$L$', '$B$'}, 'interpreter', 'latex', 'fontsize', FS_leg)

subplot(1,2,2)
p0 = plot(x_int, sum(c_ref(:,2:end).^2,2),'k:+','LineWidth',LW);
hold on
p1 = plot(x_int, sum(c_hi(:,2:end).^2,2),'-','color',c1,'LineWidth',LW);
p2 = plot(x_int, sum(c_low(:,2:end).^2,2),'--','color',c2,'LineWidth',LW);
p3 = plot(x_int, sum(c_bi(:,2:end).^2,2),'-.','color',c3,'LineWidth',LW);
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$C_p$ Variance','interpreter', 'latex', 'fontsize', FS)

new_labels = linspace(0,0.4,5);
set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend([p0,p1,p2,p3],{'Ref','$H$', '$L$', '$B$'}, 'interpreter', 'latex', 'fontsize', FS_leg)

set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,'Plots/Airfoil_mean_var','epsc')
end

% I could reverse the axis with  'YDir', 'reverse' in set(gca.. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound efficacy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/Airfoil_efficacy_1');

figure
h = pcolor(N_hi_vec, r_vec, efficacy_mat);
set(h, 'EdgeColor', 'none');
axis tight
xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Approximation Rank $r$', 'interpreter', 'latex', 'fontsize', FS)
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% title('Efficacy','interpreter', 'latex', 'fontsize', FS_leg)
set(gcf, 'Position', size_1)

if save_on == 1
    saveas(gcf,'Plots/Airfoil_efficacy','epsc')
end

figure
h = pcolor(N_hi_vec, r_vec, prob_mat);
set(h, 'EdgeColor', 'none');
axis tight
xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Approximation Rank $r$', 'interpreter', 'latex', 'fontsize', FS)
c =colorbar;
c.TickLabelInterpreter = 'latex'; 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% title('Efficacy','interpreter', 'latex', 'fontsize', FS_leg)
set(gcf, 'Position', size_1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound - single N and r 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/Airfoil_bound_results_theta_est');

% Stats that are useful: 
% eff_vec measures:  bound_23, bound_25, bound_26, bound_27, bound_36
eff_36
mean(p_33)
p_35
zeta_N
%factor_vec = [eff_vec(2), eff_vec(3)/eff_vec(2), eff_vec(4)/eff_vec(3), eff_vec(1)/eff_vec(4)];
% factor_vec

% figure
% plot(x_l,theta_vec,'-x','Color',c1, 'LineWidth',LW,'MarkerSize',MS);
% axis tight
% xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
% ylabel('$\Theta$', 'interpreter', 'latex', 'fontsize', FS)
% axis tight
% % ylim([1e-8,1])
% %xlim([1,10])
% %yticks([1e-4, 1e-2,1e0])
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% % grid on
% set(gcf,'Position',size_1)
% 
% if save_on == 1
%     saveas(gcf,'Plots/Airfoil_theta','epsc')
% end

figure
p1 = semilogy(x_l, zeta_i_1,'-x','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
hold on;
p2 = semilogy(x_l, zeta_i_2,'--o','Color',c2,...
    'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('$\zeta_{n,i}$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
legend([p1,p2],{'(22)','(24)'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,'Plots/Airfoil_zeta','epsc')
end

figure
p1 = semilogy(x_l, sqrt(err_bi_mean),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
% p2 = semilogy(x_l, sqrt(bound_20),'--s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(x_l, sqrt(bound_34),'-.x','Color',c4,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% first and last points are zero so don't plot these: 
% xlimx_lx_l(2),x_l(64)])
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p3],{'Ref Average', 'Bound (34)'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,'Plots/Airfoil_bound','epsc')
end

figure
p1 = semilogy(x_l, sqrt(err_bi_mean),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
% p2 = semilogy(x_l, sqrt(bound_20),'--s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
% p3 = semilogy(x_l, sqrt(bound_34),'-.x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
ylim([min(sqrt(err_bi_mean)), max(sqrt(bound_34))])
% first and last points are zero so don't plot these: 
% xlimx_lx_l(2),x_l(64)])
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1],{'Ref Average'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')
ylim 
if save_on == 1
    saveas(gcf,'Plots/Airfoil_bound_0','epsc')
end

figure
p1 = semilogy(x_l, sqrt(err_bi_mean),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
% p2 = semilogy(x_l, sqrt(bound_20),'--s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(x_l, sqrt(bound_34),'-.x','Color',c2,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('Location on Airfoil','interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% first and last points are zero so don't plot these: 
% xlimx_lx_l(2),x_l(64)])
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
ylim
legend([p1,p3],{'Ref Average', 'Error Bound'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,'Plots/Airfoil_bound_1','epsc')
end