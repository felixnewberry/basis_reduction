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

save_on = 0; 

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
r_string(3:end) = cellstr(num2str(r', '$B(r=%-d)$')); 

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
set(gcf, 'Position', size_1)
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)

if QoI == 1
    yl = ylim; 
    ylim([1e-3, yl(2)])
    new_labels = [1e-3, 1e-2];
    set(gca,'YTick', new_labels);
end

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)
% grid on 
set(gcf, 'Position', size_2)
yl = ylim; 
ylim([yl(1),1])
new_labels = [1e-1, 1];
set(gca,'YTick', new_labels);

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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2, p3, p4],{'$L$ $T_{\mathrm{cylinder}}$','Ref $T_{\mathrm{cylinder}}$', '$L$ $T_{x=0.2}$','Ref $T_{x=0.2}$'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'eigen'),'epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QoI mean and variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat('Results/',results_name, 'mean_var'));

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
    xticklabels({'$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'})
    xlim([-pi, pi])
end
ylabel(strcat(label_name, ' Temperature Mean'),'interpreter', 'latex', 'fontsize', FS)
legend([p0,p1,p2,p3],{'Ref','$H$','$L$', '$B$'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
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
    xticklabels({'$-\pi$','$-\pi/2$','$0$','$\pi/2$','$\pi$'})
    xlim([-pi, pi])
end
ylabel(strcat(label_name, ' Temperature Variance'),'interpreter', 'latex', 'fontsize', FS)
legend([p0,p1,p2,p3],{'Ref','$H$','$L$', '$B$'},'interpreter', 'latex', 'fontsize', FS_leg)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on

% title('Variance','interpreter', 'latex', 'fontsize', FS_axis)

set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'mean_var'),'epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Mesh             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_cylinder = [0,0,1250,225]; 

figure 

xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
xlim([-0.2, 1.0]); 
ylim([-0.1, 0.1]); 
set(gcf, 'Position', size_cylinder)



hold on
I = imread('GT_mesh_2.png'); 
h = image(xlim,ylim,I); 
uistack(h,'bottom')

set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
new_labels = [-0.1, 0, 0.1];
set(gca,'YTick', new_labels);

% hold on
% I = imread('GT_mesh_2.png'); 
% h = image(xlim,ylim,I); 
% uistack(h,'bottom')

if save_on == 1
    saveas(gcf,'Plots/GT_mesh_axis','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Colorbar for realizations             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
size_colorbar = [0,0,100,450]; 
% xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS)
% xlim([-0.2, 1.0]); 
new_x = []; 
new_labels = [2.85e2,313.75,342.5,371.25,4e2];
ylim([min(new_labels), max(new_labels)]); 
set(gca,'YTick', new_labels);
set(gca,'XTick', new_x);
set(gca, 'YAxisLocation', 'right')

set(gcf, 'Position', size_colorbar)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); %box on

ax = gca; 
% ax.YAxis.Label.String = num2str(new_labels);
ax.XAxis.Visible = 'off';



hold on
I = imread('GT_colorbar.png'); 
h = image(xlim,ylim,I); 
uistack(h,'top')



if save_on == 1
    saveas(gcf,'Plots/GT_colorbar_axis','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound efficacy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Results/GT_efficacy_cyl');
% efficacy_cy = efficacy; 
% % efficacy_cy = efficacy_mid; % Until fixed
% 
% load('Results/GT_efficacy_mid');
% efficacy_mid = efficacy; 
% 
% % Set common limits 
% lim_min = min([efficacy_cy(:); efficacy_mid(:)]);
% lim_max = max([efficacy_cy(:); efficacy_mid(:)]);
% 
% figure
% subplot(1,2,1)
% h = pcolor(N_hi_vec, r_vec, efficacy_cy);
% set(h, 'EdgeColor', 'none');
% axis tight
% xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Approximation Rank $r$', 'interpreter', 'latex', 'fontsize', FS)
% c =colorbar;
% c.TickLabelInterpreter = 'latex'; 
% caxis([lim_min, lim_max]); 
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% title('Error Bound Efficacy: Cylinder Surface','interpreter', 'latex', 'fontsize', FS_axis)
% 
% subplot(1,2,2)
% h = pcolor(N_hi_vec, r_vec, efficacy_mid);
% set(h, 'EdgeColor', 'none');
% axis tight
% xlabel('$n$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('Approximation Rank $r$', 'interpreter', 'latex', 'fontsize', FS)
% c =colorbar;
% c.TickLabelInterpreter = 'latex'; 
% caxis([lim_min, lim_max]); 
% set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% title('Error Bound Efficacy: Vertical Line','interpreter', 'latex', 'fontsize', FS_axis)
% 
% 
% 
% 
% set(gcf, 'Position', size_2)
% 
% if save_on == 1
%     saveas(gcf,'Plots/GT_efficacy','epsc')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound efficacy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat('Results/',results_name,'efficacy'));

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
    saveas(gcf,strcat('Plots/',results_name,'efficacy'),'epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound - single N and r 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat('Results/',results_name,'bound_results'));

% Stats that are useful: 
% efficacy_vec compares eqn 29, 27_sum, 42, mean(p_39) and p41
efficacy_vec
%rho_vec compares r*rho_k*2, vs expression in eqn prior - sum U_hat^2/Y...
rho_vec

figure
plot(x_h,theta_vec,'-x','Color',c1, 'LineWidth',LW,'MarkerSize',MS);
axis tight
if QoI == 0
    xlabel('$y$','interpreter', 'latex', 'fontsize', FS)
else
    xlabel('$\theta$','interpreter', 'latex', 'fontsize', FS)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'})
    xlim([-pi, pi])
end
ylabel('$\Theta$','interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'theta'),'epsc')
end

figure
plot(x_h, Y_Nh_vec,'-x','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
axis tight
if QoI == 0
    xlabel('$y$','interpreter', 'latex', 'fontsize', FS)
else
    xlabel('$\theta$','interpreter', 'latex', 'fontsize', FS)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'})
    xlim([-pi, pi])
end
ylabel('$Y$','interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
% xlim([1,10])
% yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'Y'),'epsc')
end

figure
p1 = semilogy(x_h, sqrt(err_bi_mat),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(x_h, sqrt(bound_27),'--s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(x_h, sqrt(bound_40),'-.x','Color',c4,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
if QoI == 0
    xlabel('$y$','interpreter', 'latex', 'fontsize', FS)
else
    xlabel('$\theta$','interpreter', 'latex', 'fontsize', FS)
    xticks([-pi -pi/2 0 pi/2 pi])
    xticklabels({'$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'})
    xlim([-pi, pi])
end

if QoI == 0
    yl = ylim; 
    ylim([yl(1), 1e-3])
    new_labels = [1e-4, 1e-3];
    set(gca,'YTick', new_labels);
end

ylabel('Error','interpreter', 'latex', 'fontsize', FS)
% axis tight
% first and last points are zero so don't plot these: 
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

legend([p1,p2,p3],{'True Mean','Bound (23)', 'Bound (34)'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')

if save_on == 1
    saveas(gcf,strcat('Plots/',results_name,'bound'),'epsc')
end