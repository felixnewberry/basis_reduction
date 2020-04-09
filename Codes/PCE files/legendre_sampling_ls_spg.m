

clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LW = 2;         % Line width
FS_leg = 16;    % Font size legend
FS = 28;        % Font size axis
FS_axis = 18; 
LW_axis = 2; 
size_1 = [0, 0, 410, 300]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem Switches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem = 5; 

% 0 - UB_LB_Visc
% 1 - UB_mean, UB_amplitude and UB_freq - broad scope
% 2 - UB_mean, UB_amplitude and UB_freq - narrow scope (2 versions: see problem 4 )
% 3 - UB and LB mean, amplitude and freq
% 4 - UB_mean, UB_amplitude and UB_freq - narrow scope, slightly wider
% range. 
% 5 - UB and LB mean, amplitude and freq -narrow 500
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PCE Sensitivity settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if problem==0
    % Dimension (d) and total order (p) of PC
    d = 3; % 3 for UB_LB_visc
    p = 10; % 4 for ls, 13 for spgl1

    % UB_mean UB_amp UB_freq
    N_total = 100; %140; % 140 for all bar visc
    N = 100; %140; % 80 for p=4 ls. 140, p = 13 for spgl1
    nsim_v = N_total - N;

    % UB_LB_visc
    % N_total = 100;
    % N = 70;
    % nsim_v = N_total - N;
elseif problem ==1
    % Dimension (d) and total order (p) of PC
    d = 3; 
    p = 10; % 4 for ls, 13 for spgl1

    % UB_mean UB_amp UB_freq
    N_total = 140;
    N = 140; % 80 for p=4 ls. 140, p = 13 for spgl1
    nsim_v = N_total - N;
    
elseif problem==2
    % Dimension (d) and total order (p) of PC
    d = 3; 
    p = 10; % 4 for ls, 13 for spgl1

    % UB_mean UB_amp UB_freq
    N_total = 140;
    N = 140; % 80 for p=4 ls. 140, p = 13 for spgl1
    nsim_v = N_total - N; 
elseif problem==3
    % Dimension (d) and total order (p) of PC
    d = 6; 
    p = 4; % 4 for ls, 13 for spgl1
    
%     p = 10 yields very high error.
    % p1 = 0.2451
    % p2 error 0.2452
    % p3 error: 0.3033
    % p4 error: 0.3344 (very large)
    % p5 error: 0.3175 
    % p6 error: 0.34527 (P=924)
    % p... error... 
    % UB_mean UB_amp UB_freq
    N_total = 300; % perhaps less than 300? 
    N = 300; % 80 for p=4 ls. 140, p = 13 for spgl1
    nsim_v = N_total - N; 
    
elseif problem==4
    % Dimension (d) and total order (p) of PC
    
    % really bad error... until I remove certain columns with poor time
    % convergence. 
%     What is wrong with these columns? 
    
    d = 3; 
    p = 10; % 10 - 1.78% error. 
    
    % error 17.86 % with  p = 7, 16.81
    
    N_total = 140; % perhaps less than 300? 
    N = 140; % 80 for p=4 ls. 140, p = 13 for spgl1
    nsim_v = N_total - N; 
elseif problem==5
    % Dimension (d) and total order (p) of PC
    
    % really bad error... until I remove certain columns with poor time
    % convergence. 
%     What is wrong with these columns? 
    
    d = 6; 
    p = 9; % 10 - 1.78% error. 
    
    % error 17.86 % with  p = 7, 16.81
    
    N_total = 500; % perhaps less than 300? 
    N = 500; % 80 for p=4 ls. 140, p = 13 for spgl1
    nsim_v = N_total - N; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if problem==0
    % UB_vel, LB_vel, Visc
    % QoI
    assembledData = load('assembledData_UB_LB_visc');
    
    % Random Variables
    load uniform_data_UB_LB_visc.mat;
    
    % Rename as u and y 
    u_samp = assembledData(1:N_total,1); % Use area rather than length of 
%     recirculating region. 
    y_samp = uniform_data(1:N_total,:); % samples of the input variables

elseif problem==1
    
    % assembledData_1 = load('assembledData_UBmean_UBamp_UBfreq');
    % assembledData_1 = load('assembledData_UBmean_UBamp_UBfreq_2');
    % assembledData_2 = load('assembledData_UBmean_UBamp_UBfreq_3');

    % assembledData_1 = load('assembledData_140');    
    % seperated into upper and lower area and length.
    assembledData_2 = load('assembledDataUpperLower_15_140');
    
    % load uniform_data_UBmean_UBamp_UBfreq.mat;
    % load uniform_data_UBmean_UBamp_UBfreq_2.mat;
    load uniform_data_140.mat; % I think this is 1-5. 
    
    % Upper sqrt(area)
    % u_samp = sqrt(assembledData_2(1:N_total,1)); % Use area rather than ...
%     length of recirculating region. 

    % Lower sqrt(area)
    % u_samp = sqrt(assembledData_2(1:N_total,3));

    % Both together
    u_samp = sqrt(assembledData_2(1:N_total,1) + assembledData_2(1:N_total,3));

    % u_samp = assembledData_1(1:N_total); % Use area rather than length of 
%     recirculating region. 

    % y_samp should only use UBmean, amp and freq. Identify these as columns 

    % data_array = [LBvelMax, UBvelMax, LBvelMin, UBvelMin,  LBperiod, UBperiod, ...
    %      LBfull, UBfull, LBrise, UBrise, LBfall, UBfall, Viscosity]; 

    % Uniform data 1:7: [UB_nom, UB_amp, LB_nom, LB_amp, UB_period, LB_period,
    % I should have saved UB_nom, Amp and freq

    % use for most. 
    y_samp = uniform_data(1:N_total,[1,2,5]); % samples of the input variables
    
elseif problem==2
    
    load uniform_data_UBmeanAmpFreq_34_140.mat;
    assembledData_2 = load('assembledDataUpperLower_34_140');
    u_samp = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));
%     u_samp = sqrt(assembledData_2(1:N_total,1));
    
    u_samp_upper = sqrt(assembledData_2(1:N_total,1));
    u_samp_lower = sqrt(assembledData_2(1:N_total,3));
    u_samp_comb = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));
    
    y_samp = uniform_data(1:N_total,[1,2,5]); % samples of the input variables

    1;
    
    figure(3)
    hist(u_samp_upper,30)
    xlabel('Upper Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    
    figure(4)
    hist(u_samp_lower,30)
    xlabel('Lower Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    
    figure(5)
    hist(u_samp_comb,30)
    xlabel('Combined Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on  
    
    
elseif problem==3
    
    assembledData_2 = load('assembledDataUpperLower_test5');
    assembledData_3 = load('assembledDataUpperLower_300_3');
    
    load uniform_300
%     u_samp = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));
    u_samp1 = assembledData_2;
    u_samp2 = assembledData_3;

    y_samp = uniform_data(1:N_total,[1,2,5,3,4,6]); % samples of the input variables
    % 1-6 Corresponds to UB_nom_vel, UB_nom_amp, LB_nom_vel, LB_nom_amp, UB
    % period, LB period. 
    
%     Yes, results with 0.01586 are the base case numbers. 
% Test 260 and 270... are these really zero? They did not yield files for
% the paraview errors which suggests so. 



    % Some of the samples didn't run correctly: exclude these (and row 4)

    % Uniform data 1:6: [UB_nom, UB_amp, LB_nom, LB_amp, UB_period,
    % first, third, perhaps 4th. 
    
    % LB_period, viscosity
    % I should have saved UB_nom, Amp and freq
%     [i_keep,j,s] = find(u_samp);
    
    u_test1 = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));
    u_test2 = sqrt(assembledData_3(1:N_total,1)+assembledData_2(1:N_total,3));

%     [i_remove,j,s] = find(u_samp==0);
    [i_keep1,j,s1] = find(u_test1~=u_test1(4));
    [i_keep2,j,s2] = find(u_test2~=u_test2(4));

    
%     [i_remove,j,s] = find(u_samp==u_samp(4));
%     i_remove are zeros... They have no recirc printed... check this. 

% also remove row 4. 
%     i_keep(4)=[];
    
    % Adjust vectors
    u_samp1 = u_samp1(i_keep1,:);
    u_samp2 = u_samp2(i_keep2,:);
    
    
    % look at just upper sqrt(area), or just lower area
    y_samp = y_samp(i_keep1,:); 
    
    save('u_samp1','u_samp1')
    save('u_samp2','u_samp2')

    save('y_samp','y_samp')
    
    u_samp11 = sqrt(u_samp1(:,1:2:3));
    u_samp22 = sqrt(u_samp2(:,1:2:3));
    
%     u_samp = sqrt(u_samp2(:,3));
    u_samp = sqrt(u_samp2(:,1)+u_samp2(:,3));
    
%     u_samp = sqrt(u_samp2(:,1));
    
%     u_samp_lower = sqrt(u_samp2(:,3));
    
1;

    
    
    N_total = length(u_samp);
    N=length(u_samp);

%     N_total = 282;
%     N = 250;
%     nsim_v = N_total - N;
%     

% figure
% hist(u_samp,30)
% xlabel('Recirculation Length')
% ylabel('Frequency')
% set(gca,'Fontsize', FS_axis)
elseif problem==4
    
    load uniform_24_140.mat;
    assembledData_2 = load('assembledDataUpperLower_24_140');
    % assembledData is 12 columns. correspond to upper area, upper length,
    % lower area, lower length at times t=1500, t=1000, t=500. 
    
    
    y_samp = uniform_data(1:N_total,[1,2,5]); % samples of the input variables

   
    
    y_suspect=y_samp([28,39,60,75,82,88,96],:); 
    
    results_suspect = sqrt(assembledData_2([28,39,60,75,82,88,96],[1,3,5,7,9,11]));
    % can use to look at time averaging convergence. 
    
    % Use upper and lower area together
    u_samp = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));

    % use just upper area
%     u_samp = sqrt(assembledData_2(1:N_total,1));
    
%     y_samp = uniform_data(1:N_total,[1,2,5]); % samples of the input variables
    
%     y_samp = uniform_data(1:N_total,[1,2,5,3,4,6]); % samples of the input variables
    % 1-6 Corresponds to UB_nom_vel, UB_nom_amp, LB_nom_vel, LB_nom_amp, UB
    % period, LB period. 
    
%     check params: tsize, tsteps, LBvelMax, UBvelMax, LBvelMin, UBvelMin,
%     LBperiod, UBperiod, LBfull, UBfull, LBrise, UBrise, LBfall, UBfall,
%     Viscosity. 
%     load params_24_140.mat;
% just car about columns 4, 5 8
% 4 approx 2.5-6.5, 5 0-0.9
% figure
% hist(data_array(:,5),20)

% when I look at 3-4 data: don't have readily available. compuar u_samp
% upper

% Compare convergence of upper and lower areas at time 1500, 1000, 500. 

% first: eliminate poor data: 
    assembledData_2( [28,39,60,75,82,88,96],:)=[];
    
    figure(10)
    subplot(1,2,1)
    hold on
    h1=plot(sqrt(assembledData_2(1:end,1)),'ro');
    h2=plot(sqrt(assembledData_2(1:end,5)),'kx');
    h3=plot(sqrt(assembledData_2(1:end,9)),'b*');
    hold off
    ylim([0,0.1])
    xlim([1,140])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    
    subplot(1,2,2)
    hold on
    h1=plot(sqrt(assembledData_2(1:end,3)),'ro');
    h2=plot(sqrt(assembledData_2(1:end,7)),'kx');
    h3=plot(sqrt(assembledData_2(1:end,11)),'b*');
    hold off
    ylim([0,0.1])
    xlim([1,140])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    h = legend([h3 h2 h1],{'$t_{500}$','$t_{1000}$','$t_{1500}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
% %     figure(11)
% %     subplot(1,2,1)
% %     hold on
% %     plot(sqrt(assembledData_2(1:N_total,1))-sqrt(assembledData_2(1:N_total,5)),'ro')
% %     plot(sqrt(assembledData_2(1:N_total,1))-sqrt(assembledData_2(1:N_total,9)),'kx')
% %     hold off
% %     subplot(1,2,2)
% %     hold on
% %     plot(sqrt(assembledData_2(1:N_total,3))-sqrt(assembledData_2(1:N_total,7)),'ro')
% %     plot(sqrt(assembledData_2(1:N_total,3))-sqrt(assembledData_2(1:N_total,11)),'kx')
% %     hold off
% %     
% %     % try relative error, ignore zeros. 
    
%     [i_keep3,j,s1] = find(assembledData_2(:,1) ~= 0);
%     assembledData_2 = assembledData_2(i_keep3,:);
    
    assembledData_lower = assembledData_2;

    error_upper_1000 = abs(sqrt(assembledData_2(1:end,1))...
        -sqrt(assembledData_2(1:end,5)))./sqrt(assembledData_2(1:end,1))*100;
    error_upper_500 = abs(sqrt(assembledData_2(1:end,1))...
        -sqrt(assembledData_2(1:end,9)))./sqrt(assembledData_2(1:end,1))*100;
    error_lower_1000 = abs(sqrt(assembledData_lower(:,3))...
        -sqrt(assembledData_lower(:,7)))./sqrt(assembledData_lower(:,3))*100;
    error_lower_500 = abs(sqrt(assembledData_lower(:,3))...
        -sqrt(assembledData_lower(:,11)))./sqrt(assembledData_lower(:,3))*100;
    
    figure(11)
    subplot(1,2,1)
    hold on
    h1=plot(error_upper_1000,'ro');
    h2=plot(error_upper_500,'kx');
    hold off
    ylim([0,30])
    xlim([1,140])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
%     [i_keep3,j,s1] = find(assembledData_2(:,3) ~= 0);
%     assembledData_lower = assembledData_2(i_keep3,:);
    
    subplot(1,2,2)
    hold on
    h1=plot(error_lower_1000,'ro');
    h2=plot(error_lower_500,'kx');
    hold off
    ylim([0,30])
    xlim([1,140])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    h = legend([h2 h1],{'$t_{500}$','$t_{1000}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    
    % sort by magnitude. 
    [A,sortUpper] = sort(assembledData_2(1:end,1),'descend');
    [A,sortLower] = sort(assembledData_2(1:end,3),'descend');
    
    
    figure(12)
    subplot(1,2,1)
    hold on
    h1=plot(error_upper_1000(sortUpper),'ro');
    h2=plot(error_upper_500(sortUpper),'kx');
    hold off
    ylim([0,30])
    xlim([1,140])
    xlabel('Ordered Runs', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
%     [i_keep3,j,s1] = find(assembledData_2(:,3) ~= 0);
%     assembledData_lower = assembledData_2(i_keep3,:);
    assembledData_lower = assembledData_2;
    
    subplot(1,2,2)
    hold on
    h1=plot(error_lower_1000(sortLower),'ro');
    h2=plot(error_lower_500(sortLower),'kx');
    hold off
    ylim([0,30])
    xlim([1,140])
    xlabel('Ordered Runs', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    h = legend([h2 h1],{'$t_{500}$','$t_{1000}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    
    figure(13)
    subplot(1,2,1)
    hold on
    h1=plot(sqrt(assembledData_2(1:end,1)), error_upper_1000,'ro');
    h2=plot(sqrt(assembledData_2(1:end,1)), error_upper_500,'kx');
    hold off
    ylim([0,35])
%     xlim([1,140])
    xlabel('Upper Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
%     [i_keep3,j,s1] = find(assembledData_2(:,3) ~= 0);
%     assembledData_lower = assembledData_2(i_keep3,:);
    assembledData_lower = assembledData_2;
    
    subplot(1,2,2)
    hold on
    h1=plot(sqrt(assembledData_2(1:end,3)),error_lower_1000,'ro');
    h2=plot(sqrt(assembledData_2(1:end,3)),error_lower_500,'kx');
    hold off
    ylim([0,35])
%     xlim([1,140])
    xlabel('Lower Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    h = legend([h2 h1],{'$t_{500}$','$t_{1000}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    
%     Why not look at histogram of one, then next, a way to show relative
% errors..? sort by magnitude :) make everything absolute value. too. 
% label the ybar 500 and 1000 compared to 1500

    % observations, on avergage convergence looks good, but approximately 
    
    % runs 28, 39 60, 75, 82, 88, 96 (96 and 92 have two at zero, one not... 
    
    % if I don't use these runs: 
%     adjust u_samp and y_samp


1;

    y_suspect=y_samp([28,39,60,75,82,88,96],:); 
    % doesn't seem to be anything odd... ho hum. 
    
    u_samp_comb = sqrt(assembledData_2(1:end,1)+ assembledData_2(1:end,3));
    u_samp_upper = sqrt(assembledData_2(1:end,1));
    u_samp_lower = sqrt(assembledData_2(1:end,3));
    
    u_samp( [28,39,60,75,82,88,96],:)=[];
    y_samp( [28,39,60,75,82,88,96],:)=[];
    
    N_total = length(u_samp);
    N=length(u_samp);
    
    figure(3)
    hist(u_samp_upper,30)
    xlabel('Upper Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    
    figure(4)
    hist(u_samp_lower,30)
    xlabel('Lower Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    
    figure(5)
    hist(u_samp_comb,30)
    xlabel('Combined Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on  
    
elseif problem==5
    
    load uniform_500.mat;
    
    % load area
% %     assembledData_2 = load('assembledDataUpperLower_500');
% %     % assembledData is 12 columns. correspond to upper area, upper length,
% %     % lower area, lower length at times t=2000, t=1500, t=1000. 
% %     
% %     % Use upper and lower area together
% % %    u_samp = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));
% % 
% %     % use just upper area
% % %     u_samp = sqrt(assembledData_2(1:N_total,1));
% %     
% % % use just lower area. 
% %     u_samp = sqrt(assembledData_2(1:N_total,1));
% % 
% %     u_samp_upper = sqrt(assembledData_2(1:N_total,1));
% %     u_samp_lower = sqrt(assembledData_2(1:N_total,3));
% %     u_samp_comb = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,3));
% %     
% %     u_samp_all = [u_samp_upper, u_samp_lower, u_samp_comb]; 

    % load area wrt dwall. 
    assembledData_2 = load('assembledDwall_500');
    % assembledData is 6 columns. correspond to dwall area, upper, at times t=2000, t=15000, t=1000 followed by lower 
    
    % Use upper and lower area together
   u_samp = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,4));

    % use just upper area
%     u_samp = sqrt(assembledData_2(1:N_total,1));
    
% use just lower area. 
%     u_samp = sqrt(assembledData_2(1:N_total,4));

    u_samp_upper = sqrt(assembledData_2(1:N_total,1));
    u_samp_lower = sqrt(assembledData_2(1:N_total,4));
    u_samp_comb = sqrt(assembledData_2(1:N_total,1)+assembledData_2(1:N_total,4));
    
    u_samp_all = [u_samp_upper, u_samp_lower, u_samp_comb]; 
    
    
    % convert to compatible form for convergence plots
    temp_save = assembledData_2;
    assembledData_2 = zeros(N_total, 12); 
    
    assembledData_2(:,[1,5,9,3,7,11])= temp_save;
    
    
    y_samp = uniform_data(1:N_total,[1,2,5,3,4,6]); % samples of the input variables
    
    1;
    % Find zeros to show on contour plot
    [i_zeros,j,s1] = find(u_samp_comb(:,1) == 0);
    
    y_zeros= y_samp(i_zeros,:);
    u_samp_all;
    %     save('u_samp_all','u_samp_all')

%     save('y_samp','y_samp')
%     save('y_zeros','y_zeros')
    % 1-6 Corresponds to UB_nom_vel, UB_nom_amp, LB_nom_vel, LB_nom_amp, UB
    % period, LB period. 
    
%     check params: tsize, tsteps, LBvelMax, UBvelMax, LBvelMin, UBvelMin,
%     LBperiod, UBperiod, LBfull, UBfull, LBrise, UBrise, LBfall, UBfall,
%     Viscosity. 
%     load params_24_140.mat;
% just car about columns 4, 5 8
% 4 approx 2.5-6.5, 5 0-0.9
% figure
% hist(data_array(:,5),20)

% when I look at 3-4 data: don't have readily available. compuar u_samp
% upper

% Compare convergence of upper and lower areas at time 1500, 1000, 500. 

% first: eliminate poor data: 
%     assembledData_2( [28,39,60,75,82,88,96],:)=[];
    
    figure(10)
    subplot(1,2,1)
    hold on
    h1=plot(sqrt(assembledData_2(1:N_total,1)),'ro');
    h2=plot(sqrt(assembledData_2(1:N_total,5)),'kx');
    h3=plot(sqrt(assembledData_2(1:N_total,9)),'b*');
    hold off
    ylim([0,0.1])
    xlim([1,500])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    
    subplot(1,2,2)
    hold on
    h1=plot(sqrt(assembledData_2(1:N_total,3)),'ro');
    h2=plot(sqrt(assembledData_2(1:N_total,7)),'kx');
    h3=plot(sqrt(assembledData_2(1:N_total,11)),'b*');
    hold off
    ylim([0,0.1])
    xlim([1,500])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    h = legend([h3 h2 h1],{'$t_{1000}$','$t_{1500}$','$t_{2000}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    assembledData_lower = assembledData_2;
    error_upper_1000 = abs(sqrt(assembledData_2(1:N_total,1))...
        -sqrt(assembledData_2(1:N_total,5)))./sqrt(assembledData_2(1:N_total,1))*100;
    error_upper_500 = abs(sqrt(assembledData_2(1:N_total,1))...
        -sqrt(assembledData_2(1:N_total,9)))./sqrt(assembledData_2(1:N_total,1))*100;
    error_lower_1000 = abs(sqrt(assembledData_lower(:,3))...
        -sqrt(assembledData_lower(:,7)))./sqrt(assembledData_lower(:,3))*100;
    error_lower_500 = abs(sqrt(assembledData_lower(:,3))...
        -sqrt(assembledData_lower(:,11)))./sqrt(assembledData_lower(:,3))*100;
    
    figure(11)
    subplot(1,2,1)
    hold on
    h1=plot(error_upper_1000,'ro');
    h2=plot(error_upper_500,'kx');
    hold off
    ylim([0,10])
    xlim([1,500])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{2000}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
%     [i_keep3,j,s1] = find(assembledData_2(:,3) ~= 0);
%     assembledData_lower = assembledData_2(i_keep3,:);
    assembledData_lower = assembledData_2;
    
    subplot(1,2,2)
    hold on
    h1=plot(error_lower_1000,'ro');
    h2=plot(error_lower_500,'kx');
    hold off
    ylim([0,10])
%     xlim([1,500])
    xlabel('Run', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{1500}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    h = legend([h2 h1],{'$t_{1000}$','$t_{1500}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    
    % sort by magnitude. 
    [A,sortUpper] = sort(assembledData_2(1:N_total,1),'descend');
    [A,sortLower] = sort(assembledData_2(1:N_total,3),'descend');
    
    
    figure(12)
    subplot(1,2,1)
    hold on
    h1=plot(error_upper_1000(sortUpper),'ro');
    h2=plot(error_upper_500(sortUpper),'kx');
    hold off
    ylim([0,10])
%     xlim([1,500])
    xlabel('Ordered Runs', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{2000}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
%     [i_keep3,j,s1] = find(assembledData_2(:,3) ~= 0);
%     assembledData_lower = assembledData_2(i_keep3,:);
    assembledData_lower = assembledData_2;
    
    subplot(1,2,2)
    hold on
    h1=plot(error_lower_1000(sortLower),'ro');
    h2=plot(error_lower_500(sortLower),'kx');
    hold off
    ylim([0,10])
    xlim([1,500])
    xlabel('Ordered Runs', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Relative Error to $t_{2000}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    h = legend([h2 h1],{'$t_{1000}$','$t_{1500}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    
    figure(13)
    subplot(1,2,1)
    hold on
    h1=plot(sqrt(assembledData_2(1:N_total,1)), error_upper_1000,'ro');
    h2=plot(sqrt(assembledData_2(1:N_total,1)), error_upper_500,'kx');
    hold off
    ylim([0,10])
%     xlim([1,140])
%     xlabel('Upper Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    xlabel('Upper Integrated dwall [$m^{-3}$] ', 'interpreter', 'latex', 'fontsize', FS)

    ylabel('Relative Error to $t_{2000}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
%     [i_keep3,j,s1] = find(assembledData_2(:,3) ~= 0);
%     assembledData_lower = assembledData_2(i_keep3,:);
    assembledData_lower = assembledData_2;
    
    subplot(1,2,2)
    hold on
    h1=plot(sqrt(assembledData_2(1:N_total,3)),error_lower_1000,'ro');
    h2=plot(sqrt(assembledData_2(1:N_total,3)),error_lower_500,'kx');
    hold off
    ylim([0,10])
%     xlim([1,140])
%     xlabel('Lower Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    xlabel('Lower Integrated dwall [$m^{-3}$] ', 'interpreter', 'latex', 'fontsize', FS)

    ylabel('Relative Error to $t_{2000}$', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    h = legend([h2 h1],{'$t_{1000}$','$t_{1500}$'});
    set(h, 'interpreter', 'latex', 'fontsize', FS_leg)
    
    
%     Why not look at histogram of one, then next, a way to show relative
% errors..? sort by magnitude :) make everything absolute value. too. 
% label the ybar 500 and 1000 compared to 1500

    % observations, on avergage convergence looks good, but approximately 
    
    % runs 28, 39 60, 75, 82, 88, 96 (96 and 92 have two at zero, one not... 
    
    % if I don't use these runs: 
%     adjust u_samp and y_samp
    
    figure(3)
    hist(u_samp_upper,30)
%     xlabel('Upper Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    xlabel('Upper Integrated dwall [$m^{-3}$] ', 'interpreter', 'latex', 'fontsize', FS)

    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on
    
    figure(4)
    hist(u_samp_lower,30)
    xlabel('Lower Integrated dwall [$m^{-3}$] ', 'interpreter', 'latex', 'fontsize', FS)

%     xlabel('Lower Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on    
    
    figure(5)
    hist(u_samp_comb,30)
    xlabel('Combined Integrated dwall [$m^{-3}$] ', 'interpreter', 'latex', 'fontsize', FS)

%     xlabel('Combined Recirculation Length', 'interpreter', 'latex', 'fontsize', FS)
    ylabel('Frequency', 'interpreter', 'latex', 'fontsize', FS)
    set(gca,'Fontsize', FS_axis)
    grid on  
    
end

%%% Look at error of upper and lower areas... 

% load uniform_data_UBmeanAmpFreq_34_140.mat;

% narrow scope

% rel_error = abs((assembledData_1 - assembledData_2)./assembledData_2*100);
% figure(1)
% histogram(rel_error,20)

% data_cols = [assembledData_1,assembledData_2,assembledData_1 - assembledData_2];
% [a,b] = max(rel_error);

% yeilds approximatley 8 % error. Quite high. Rel error doing an extra 1000
% samples has bout 20% of samples in range 0.5-2%. 
% RMS is negligable for lower seperation, and most of upper seperation.
% About 1-1.4 % for limited area of upper circulation. 

%take average y_bar? 
% assembledData_1 = mean((assembledData_2 + assembledData_1),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sensitivity Calculation - make neater. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tolerance on residual used in spgl1 solver
% sigma = 1e-5;

% Important parameters for spgl1 toolbox
% verbosity = 0 suppresses spgl1 command line output
% opts = spgSetParms('iterations',10000,'verbosity',0);
opts = spgSetParms('iterations',10000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);

% Reference PCE coefficents via l_1 minimization
% spg_mmv is a function in the spgl1 toolbox -Solve multi-measurement basis pursuit denoise
% c_ref = spg_mmv(psi_ref,u_ref,sigma*norm(u_ref),opts);

% Load the PCE idex matrix
index_pc = nD_polynomial_array(d,p);

% size of the basis
P = size(index_pc,1);

% Solve the l1 minimization problem min || c ||_1 s.t. ||Psi c - u_samp||_2^2 < delta 

% Form the measurement matrix Ps
for isim=1:N
    crow = piset(y_samp(isim,:),index_pc);
    Psi(isim,:) = crow;
end

% Solve for the solution using least squares

% c_ls = Psi\u_samp(1:N,1);

% or solve with spgl1:
% Reference PCE coefficents via l_1 minimization
% spg_mmv is a function in the spgl1 toolbox -Solve multi-measurement basis pursuit denoise
% c_ref = spg_bpdn(psi_ref,u_ref,sigma,opts);

% % Psi - dictionary of polynomials
% % 




weights = get_matrix_weights(Psi);
Psiweights = Psi*weights;
% 
% sigma is truncation error of PCE (approximated)
sigma =  cross_val_sigma(Psi,u_samp(1:N,1));
c_spg = weights*spg_bpdn(Psiweights,u_samp(1:N,1),sigma*norm(u_samp(1:N,1)),opts);

% % c_spg = SP(Psi,0,u_samp,0,'random');
% % val_error = norm(Psi*c_spg - u_samp)/norm(u_samp)

% % c_spg = OMP(Psi,-Inf,u_samp,0);

val_error = norm(Psi*c_spg - u_samp)/norm(u_samp)

% c_spg = c_ls;
% % save('c_spg_UB_LB_500_comb_dwall','c_spg')
% % save('c_spg_UB_133_comb','c_spg')

% Plot the coefficients 
figure(2)
loglog(abs(c_spg),'ro'); 


% legend('l2 minimization')

% Relative mean square error
% err_rms = norm(c_spg-c_ref(1:P))/norm(c_ref(1:P)); disp(err_rms);


%% Start validation

% clear Psi;

% for isim=1:nsim_v
%     crow = piset(y_samp(isim+N,:),index_pc);
%     Psi(isim,:) = crow(1:P);
% end
% error_val = norm(u_samp(N+1:end,1)-Psi*c_spg)/norm(u_samp(N+1:end,1))

%% Comput the Sobol' indices

% Call the Sobol index code
s = zeros(1,d);

for j=1:d    
    s(j) = sobol_pce( c_spg, index_pc(:,j));
end

% Plot Sobol' indices - make generic for number of variables... 

% number of input variables is d: 

labels=cell(1,d);
labels_blank=cell(1,d);
for i = 1:d
    labels{i}=strcat('\xi_', num2str(i));
    labels_blank{i}=strcat('');
end
% labels_1 = sym('xi_',[d 1]);

% set some labels to zero, ie 3 and 6
% labels{3}='';
% labels{6}='';

% A_cell = cellstr(str2num(labels_1));

% For Problem 3:
% Error using heaviside (line 22)
% Invalid data type. Argument must be single, double, or sym.
% s is 274*1... all ones... doesn't really seem right. 

% make a color bar on right hand side of plot.  

explode = heaviside(s);
figure(4)
h = pie(s,explode,labels_blank);

% obtain colormap of figure
cmap = colormap;

colormap(parula(d));

c = colorbar; 
tick_position= linspace(1,d,2*d+1);
% tick_position= linspace(1,d-1,2*d+1);

tick_position = tick_position(2:2:end); 
%tick_position=tick_position(1:d)+(tick_position(2)-tick_position(1))/2;
c.Ticks = tick_position;
c.TickLabels = labels;
w = c.LineWidth';
c.LineWidth = 1.5; 
c.Box = 'on';
% change width and size of colorbar:
x1=get(gca,'position');
x=get(c,'Position');

bar_width = 0.05; 
% bar_height = 
x(3)=bar_width; % colorbar width

% y_mid = 0.5*(x(2)+x(4));
% x(2)= y_mid-bar_width*d;
x(4)= 2*bar_width*d; % colorbar height? 
x(2) = (x1(2)+x1(4))/2-bar_width*d;
set(c,'Position',x);
set(gca,'position',x1);


set(gca, 'FontSize',FS);
% c.Ticks.FontSize=FS;
% c.TickLabelInterpreter ='latex','\fontsize{FS} labels';

% d colours... 
% colormap([255 244 102; %// yellow
%         112 193 129;      %// green
%         36 123 160]./255) %// blue

% % Set up label sizes
% for i_text = 1:d
%     h(2*i_text).FontSize=18;
% end
% 
% set(gcf, 'Position', size_1)


% set(gca,'Fontsize', FS_axis);
% xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
% ylabel('$\lambda_i$', 'interpreter', 'latex', 'fontsize', FS)
% set(gca,'linewidth',LW_axis);
% box on
% grid on

% Plot Sobol' indices - make generic for number of variables... 
% labels = {' \xi_1 ','    \xi_2 ',' \xi_3 ',' \xi_4 ','    \xi_5 ',' \xi_6    '};
% explode = heaviside(s);
% figure(3)

% plot first order sobal index. 

% % p = pie(s,explode,labels);
% % colormap([255 244 102; %// yellow
% %         112 193 129;      %// green
% %         36 123 160; %// blue
% %         255 244 102; %// yellow
% %         112 193 129;      %// green
% %         36 123 160]./255) %// blue
% % 
% %     t1 = p(2); 
% % t2 = p(4); 
% % t3 = p(6); 
% % 
% % t1.FontSize = 18;
% % t2.FontSize = 18; 
% % t3.FontSize = 18; 

% % % set(gca,'Fontsize', FS_axis);
% % % xlabel('index $i$', 'interpreter', 'latex', 'fontsize', FS)
% % % ylabel('$\lambda_i$', 'interpreter', 'latex', 'fontsize', FS)
% % % set(gca,'linewidth',LW_axis);
% % % box on
% % set(gcf, 'Position', size_1)
% % % grid on

% figure(2)
% plot_n = 4;
% subplot(1,2,1)
% hist(data_array_300(:,plot_n),30) 
% subplot(1,2,2)
% hist(data_array_200(:,plot_n),30,'facecolor','r','facealpha',0.5) 
% 
% figure(3)
% plot_n = 6;
% subplot(1,2,1)
% hist(data_array_300(:,plot_n),30) 
% subplot(1,2,2)
% hist(data_array_200(:,plot_n),30,'facecolor','r','facealpha',0.5) 

% data_array = [ Tsize, Tsteps, LBvelMax, UBvelMax, LBvelMin, UBvelMin,  LBperiod, UBperiod, ...
%      LBfull, UBfull, LBrise, UBrise, LBfall, UBfall, Viscosity]; 

% Inpspecting data of runs 4 and 6. Seem contrary. 

load('/Users/felixnewberry/Google Drive/1_Research/9_SSE/2D_diffuser/Global_sensitivity/Data/params_24_140.mat')
data_test_4 = data_array;

load('/Users/felixnewberry/Google Drive/1_Research/9_SSE/2D_diffuser/Global_sensitivity/Data/params_500.mat')
data_test_6 = data_array;

% Data is in max and min form. Convert to mean and amp:

t4_mean = (data_test_4(:,4)+data_test_4(:,6))/2;
t4_amp = data_test_4(:,4)-t4_mean;
t6_mean = (data_test_6(:,4)+data_test_6(:,6))/2;
t6_amp = data_test_6(:,4)-t6_mean;

% identify runs that are 0 in Upper or lower: 
assembledData_4 = load('assembledDataUpperLower_24_140');
u_samp4_upper = sqrt(assembledData_4(1:end,1));
u_samp4_lower = sqrt(assembledData_4(1:end,3));
assembledData_6 = load('assembledDataUpperLower_500');
u_samp6_upper = sqrt(assembledData_6(1:end,1));
u_samp6_lower = sqrt(assembledData_6(1:end,3));

[i_zeros4,j,s1] = find(u_samp4_upper(:,1) == 0);
[i_zeros6,j,s1] = find(u_samp6_upper(:,1) == 0);

figure
hold on
p1 = plot(t4_mean, t4_amp, 'ro'); 
p2 = plot(t6_mean, t6_amp, 'bx'); 
p3 = plot(t4_mean(i_zeros4), t4_amp(i_zeros4), 'ko'); 
p4 = plot(t6_mean(i_zeros6), t6_amp(i_zeros6), 'kx'); 
hold off

% Inspect run that has mean 3.879, amp 2.861 for just upper blower. Compare
% to 3.864, 2.861 with both blowers. 

% [i_test4,j,s1] = find(t4_mean == 3.879);
% 
% for i = 1:15
%     [inputs_97(i) inputs_353(i)]
% end
