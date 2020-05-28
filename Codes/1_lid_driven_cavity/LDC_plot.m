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

save_on = 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Load data                  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Results/LDC_results.mat')

% Vector of strings for r plots
r_symbol = {'-.+','-.*','-.s','-.d'}; 

r_string = cell(length(r)+2,1);
r_string(1:2,:) = {'$H$'; '$L$'};
r_string(3:end) = cellstr(num2str(r', 'B(r=%-d)')); 

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
xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
ylabel('$\lambda_i$', 'interpreter', 'latex', 'fontsize', FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% grid on
set(gcf,'Position',size_1)
xlim([1,5])
legend([p1,p2],{'L','Ref'},'interpreter', 'latex', 'fontsize', FS_leg,'Location','NorthEast')

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
ylabel('Average relative errors in variance','Fontsize',FS)
axis tight
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
legend(r_string,'interpreter', 'latex', 'fontsize', FS_leg)
% grid on 
set(gcf, 'Position', size_2)

if save_on == 1
    saveas(gcf,'Plots/LDC_N_r','epsc')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LDC QoI
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

xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
axis tight
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis); %box on
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
xlabel('x','interpreter', 'latex', 'fontsize', FS)
ylabel('y','interpreter', 'latex', 'fontsize', FS)
xlim([0,1]); ylim([0,1]);
new_labels = linspace(0, 1, 3);
set(gca,'XTick', new_labels); set(gca,'YTick', new_labels);
set(gca,'Fontsize', FS_axis, 'linewidth',LW_axis);box on
% set(gcf, 'Position', size_1)
set(gcf, 'Position', size_square)

pbaspect([1 1 1])
if save_on ==1
    saveas(gcf,'Plots/LDC_Geom_stream','epsc')
end

% Should I make a note of the variance?? 
% 100 repetitions