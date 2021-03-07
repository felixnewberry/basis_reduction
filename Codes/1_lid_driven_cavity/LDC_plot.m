clear all 
close all 
clc

% Plot LDC

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

load('Results/LDC_results_spg.mat')

% Vector of strings for r plots
r_symbol = {'-.+','-.*','-.s','-.d'}; 

r_string = cell(length(r)+2,1);
r_string(1:2,:) = {'$H$'; '$L$'};
r_string(3:end) = cellstr(num2str(r', '$B(r=%-d)$')); 

r_string_variance = r_string([1,3:end]); 

% r_plot = [1 2 3];
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
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
xlim([1,5])
legend([p1,p2],{'$L$','Ref'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

if save_on == 1
    saveas(gcf,'Plots/LDC_eigen','epsc')
end

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
for i_r = 1:length(r)
    semilogy(N_hi,(mean_mean_bi_err(i_r,:)),...
        r_symbol{i_r},'Color',c3, 'LineWidth',LW,'MarkerSize',MS)
end

ylabel('Average Relative Error in Mean','interpreter', 'latex', 'fontsize', FS)
xlabel('$n$','interpreter','latex','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex');box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)
% title('Mean','interpreter', 'latex', 'fontsize', FS_axis)

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
ylabel('Average Relative Error in Variance','interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex');box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)
% title('Variance','interpreter', 'latex', 'fontsize', FS_axis)
set(gcf, 'Position', size_2)


if save_on == 1
    saveas(gcf,'Plots/LDC_N_r','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LDC QoI - flow field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plot contour and then streamlines

load('LDC_data/stream_coords')
x_stream = x; 
y_stream = y; 

load 'LDC_data/u_out_plot.mat'

n_cell_x = 128; %CHANGE TO 128.
nx = n_cell_x;
x = linspace(0,1,n_cell_x+1);


x = (x-0.5)*2;
x = 0.5*(cos(pi*(x-1)/2)+1);
y = x;

u = u_out(1:end/2);
v = u_out(end/2+1:end);


u_mat = reshape(u,nx+1,nx+1)'; 
v_mat = reshape(v,nx+1,nx+1)'; 

u_mag = sqrt(u_mat.^2 +v_mat.^2);

x_1 = linspace(0,1.0,11);
y_1 = 0.5*ones(1, length(x_1)); 


step_quiv = 10; 

figure
[c,h]=contourf(x,y,u_mag);
set(h, 'edgecolor','none');
hold on

quiver(x(1:step_quiv:end),y(1:step_quiv:end),u_mat(1:step_quiv:end,1:step_quiv:end), v_mat(1:step_quiv:end,1:step_quiv:end),'Color',c3)

% plot qoi
p1 = plot(x_1, y_1,'--','color',c2,'LineWidth',LW+1);

X_arrow = [0.35 0.7];
Y_arrow = [0.95   0.95];
hh = annotation('arrow',X_arrow,Y_arrow,'Color','r');
set(hh, 'LineWidth', LW)

xlabel('$x$','interpreter', 'latex', 'fontsize', FS+5)
ylabel('$y$','interpreter', 'latex', 'fontsize', FS+5)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); %box on
set(gcf, 'Position', size_square)
pbaspect([1 1 1])

if save_on ==1
    saveas(gcf,'Plots/LDC_Geom_contour','epsc')
end

% try quiver... 
% select step size: 



% Option1, standard. 
figure
h = streamslice(x,y, u_mat, v_mat, 1, 'cubic'); %, 0.9,0.1, [0.1, 8000])
set( h, 'Color', c1 )
hold on
set(h, 'LineWidth', LW) %LW/2 ?

% plot qoi
p1 = plot(x_1, y_1,'--','color',c2,'LineWidth',LW+1);

X_arrow = [0.35 0.7];
Y_arrow = [0.95   0.95];
hh = annotation('arrow',X_arrow,Y_arrow,'Color','r');
set(hh, 'LineWidth', LW)

hold off
xlabel('$x$','interpreter', 'latex', 'fontsize', FS+5)
ylabel('$y$','interpreter', 'latex', 'fontsize', FS+5)
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); %box on
% set(gcf, 'Position', size_1)
set(gcf, 'Position', size_square)

pbaspect([1 1 1])
if save_on ==1
    saveas(gcf,'Plots/LDC_Geom_stream','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_l = 4; 
n_h = 64; 

x_l  = linspace(0,1,n_l); 
x_h  = linspace(0,1,n_h); 


% concentrate points close to sides
x_l = (x_l - 0.5).*2; 
x_l = 0.5.*(cos(pi.*(x_l-1) / 2) + 1);

[xx_l, yy_l ] = meshgrid(x_l); 

x_h = (x_h - 0.5).*2; 
x_h = 0.5.*(cos(pi.*(x_h-1) / 2) + 1);

[xx_h, yy_h ] = meshgrid(x_h); 

% plot mesh.. 
% plot
figure
plot(xx_l, yy_l, 'k','LineWidth',LW/2)
hold on; 
plot(yy_l, xx_l, 'k','LineWidth',LW/2)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS+5)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS+5)
% grid on; 
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); %box on
axis tight;

set(gcf,'Position',size_square)

if save_on ==1
    saveas(gcf,'Plots/LDC_mesh_low','epsc')
end

figure
plot(xx_h, yy_h, 'k','LineWidth',LW/4)
hold on; 
plot(yy_h, xx_h, 'k','LineWidth',LW/4)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS+5)
ylabel('$y$', 'interpreter', 'latex', 'fontsize', FS+5)
% grid on; 
set(gca,'Fontsize', FS_axis+5, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); %box on
axis tight;
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gcf,'Position',size_square)


if save_on ==1
    saveas(gcf,'Plots/LDC_mesh_high','epsc')
end

1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% QoI mean and variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/LDC_qoi_mean_var');

figure
subplot(1,2,1)
p0 = plot(x_h, c_ref(:,1),'k:+','LineWidth',LW);
hold on
p1 = plot(x_h, c_hi(:,1),'-','color',c1,'LineWidth',LW);
p2 = plot(x_h, c_low(:,1),'--','color',c2,'LineWidth',LW);
p3 = plot(x_h, c_bi(:,1),'-.','color',c3,'LineWidth',LW);
xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
ylabel('Vertical Velocity Mean','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% title('Mean','interpreter', 'latex', 'fontsize', FS_axis)
legend([p0,p1,p2,p3],{'Ref','$H$','$L$', '$B$'},'interpreter', 'latex', 'fontsize', FS_leg)

subplot(1,2,2)
p0 = plot(x_h, sum(c_ref(:,2:end).^2,2),'k:+','LineWidth',LW);
hold on
p1 = plot(x_h, sum(c_hi(:,2:end).^2,2),'-','color',c1,'LineWidth',LW);
p2 = plot(x_h, sum(c_low(:,2:end).^2,2),'--','color',c2,'LineWidth',LW);
p3 = plot(x_h, sum(c_bi(:,2:end).^2,2),'-.','color',c3,'LineWidth',LW);
xlabel('$x$','interpreter', 'latex', 'fontsize', FS)
ylabel('Vertical Velocity Variance','interpreter', 'latex', 'fontsize', FS)
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% title('Variance','interpreter', 'latex', 'fontsize', FS_axis)
legend([p0,p1,p2,p3],{'Ref','$H$','$L$', '$B$'},'interpreter', 'latex',...
    'fontsize', FS_leg, 'Location','northwest')

set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,'Plots/LDC_mean_var','epsc')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound efficacy 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/LDC_efficacy_2');
% load('LDC_bound_results_test_t');

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
    saveas(gcf,'Plots/LDC_efficacy','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bound - single N and r 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('Results/LDC_bound_results');
load('LDC_bound_results_test_t');
% Stats that are useful: 
% efficacy_vec compares eqn 29, 27_sum, 42, mean(p_39) and p41
% Translates to... 
efficacy_vec
%rho_vec compares r*rho_k*2, vs expression in eqn prior - sum U_hat^2/Y...
rho_vec

figure
plot(x_h,theta_vec,'-x','Color',c1, 'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\Theta$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

if save_on == 1
    saveas(gcf,'Plots/LDC_theta','epsc')
end

figure
semilogy(x_h, Y_Nh_vec,'-x','Color',c1,...
    'LineWidth',LW,'MarkerSize',MS);
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$Y$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)
if save_on == 1
    saveas(gcf,'Plots/LDC_Y','epsc')
end

figure
p1 = semilogy(x_h, sqrt(err_bi_mat),'-o','Color',c1,'LineWidth',LW,'MarkerSize',MS);
hold on
p2 = semilogy(x_h, sqrt(bound_16),'--s','Color',c3,'LineWidth',LW,'MarkerSize',MS);
p3 = semilogy(x_h, sqrt(bound_34),'-.x','Color',c4,'LineWidth',LW,'MarkerSize',MS);
% p4 = semilogy(x_h, sqrt(bound_24),'-.<','Color',c5,'LineWidth',LW,'MarkerSize',MS);
hold off
axis tight
xlabel('$x$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('Error', 'interpreter', 'latex', 'fontsize', FS)
axis tight
% first and last points are zero so don't plot these: 
xlim([x_l(2),x_l(64)])
% ylim([1e-8,1])
%xlim([1,10])
%yticks([1e-4, 1e-2,1e0])
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis,'TickLabelInterpreter','latex'); box on
% grid on
set(gcf,'Position',size_1)

% legend([p1,p2,p3, p4],{'Reference Average','Bound (16)', 'Bound (34)',  'Bound (24ish)'}...
%     ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')
legend([p1,p2,p3],{'Ref Average','Bound (16)', 'Bound (34)'}...
    ,'interpreter', 'latex', 'fontsize', FS_leg,'Location','SouthEast')

if save_on == 1
    saveas(gcf,'Plots/LDC_bound','epsc')
end
%



efficacy_vec
%rho_vec compares r*rho_k*2, vs expression in eqn prior - sum U_hat^2/Y...
rho_vec
