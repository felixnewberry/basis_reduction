clear all 
close all
clc

tic
% pause(0.1)
% display('test')

%%% Felix Newberry: Code BR (basis reduction for familiarity)

%This is the main script from which other scripts will be called

% %%%%%%%%%%%%%%%%%%%%% Task Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select model: only 1. 
run_LDC = 1;        % run lid driven cavity
run_airfoil = 0;    % run airfoil
run_LIB = 0;        % run Lithium Ion Battery
run_GT = 0;        % run Gas Turbine

% Trying to introduce spg as opposed to spg mmv, yet to succeed... 

%% Problem specific section
% Problem specific in the sense that inputs and context of problem will
% change. The structure of this secion should remain constant. 

% Purpose
% 1)load high and low fidelity response data and rename as u_ref (for high
%   fid data) and u_low (for low fid data). 
% 
% 2)Load UNIFORMLY distributed random variables and normalize them between [
%   -1 and 1]. Re-name as xi. 
%   Load response data for high fid model. (all available samples)
% 
% 3)Select your PCE total order, approximation rank, the range of number of
%   high fid samples. 

%%% Lid driven cavity (1d, may have to edit some things... list them
%%% clearly)

% Low fidelity data

% Contents for lid driven cavity data:
% -kinematic viscosity nu 200 x 1
% -lid velocity v 200 x 1
% -data points u n_points x 200. ie v velocity at ponts along y = 0.5 slice

% Q: Both high and low have the same pertubations. Is this necessary? 
% %%%%%%%%%%%%%%%%%%%%% Airfoil %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_airfoil ==1
    load lowFiResults_6d2.mat;

    % load response data for low fidelity model (all available samples)
    load highFiResults_6d2.mat;

    % I only want to look at the first few grid points. there are 128 total
    % available
    gridpt = [1:128];

    % u_ref should be Number of samples x spatial/temporal dimension of problem
    u_ref = highFiResults(gridpt,:)';


    % u_low should be Number of samples x spatial/temporal dimension of problem
    u_low = lowFiResults(gridpt,:)';

    % load random variables used to generate u
    load Inputs_6d2.mat;

    load('x_locations')
    x_l = x_locations; 
    x_h = x_locations;

%     normalize to be U~[-1,1] and re-label as xi.

    % xi_ref should be of dimension: Number of reference samples x stochastic
    % dimension d
    xi_ref =[(1/3)*alpha',(4*c-4)',(50*m-1)',(10*(mach-.2)-2)',(10*p -4)',((t - .1)*50 - 1)'];

    % x_low should be of dimension: number of low fidelity samples X stochastic
    % dimension d
    xi_low = xi_ref;
end
% %%%%%%%%%%%%%%%%%%%%% LDC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_LDC == 1
%     load Results_lid_8

    % Store data
    % lowFiResults = u_vec; 

    load lowFiResults_LDC

    % % Lid driven cavity specific adjustment
    % Re = 100; 
    % nu_av = 1/Re; 
    % 
    % 
    % % Range of pertubation
    % delta_max_nu = nu_av*0.05; %plus or minus delta_max on Left Hand boundary condition
    % delta_max_v = 0.1; 
    % 
    % % Store pertubations/ RVs. Size: n_samp x n_stochastic parameters
    % xi_low = [nu, v]; 
    % xi_low = [xi_low(:,1)/delta_max_nu, xi_low(:,2)/delta_max_v]; 

    load xi_low_LDC

    % High fidelity data
    % load Results_lid_32

    load highFiResults_LDC.mat

    % highFiResults = u_vec; 

    % % Store pertubations/ RVs. Size: n_samp x n_stochastic parameters
    % xi_high = [nu, v]; 
    % xi_high = [xi_high(:,1)/delta_max_nu, xi_high(:,2)/delta_max_v]; 

    load xi_high_LDC.mat

    % data points are 65 along y = 0.5

    % Load random inputs
    % Already normalized between -1 and 1
    % xi_ref=xi_high=xi_low for LDC, is this necessary for error bound? 
    xi_ref = xi_high; 


    % Grid points of interest
    gridpt_l = 1:size(lowFiResults(:,1)); 
    gridpt_h = 1:size(highFiResults(:,1)); 


    % u_ref should be number of samples x spatial/temporal dimension of the
    % problem. size: n_samples x n_gridpoints
    u_ref = highFiResults(gridpt_h,:)';

    % u_low should be number of samples x spatial/temporal dimension of the
    % problem. size: n_samples x n_gridpoints
    u_low = lowFiResults(gridpt_l,:)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Lithium Ion Battery
if run_LIB ==1
    
    % load data, make work. Nail it. Revise PCE, Revise KL, Revise basis
    % produciton. nail slides. 


    % Contents for lithium ion battery 
    % c_c has the low-fidelity data matrix that has 587 samples of concentration 
    % computed over the battery cell and on a grid of size 50
    % grid_coarse.mat contains the coordinates of of the coarse grid points
    % c_f and grid_fine are the high-fidelity counterparts of the above. 
    % rv.mat has the realizations of the random inputs (17 inputs in total).

    % load and store low fidelity results
    load c_c.mat
    lowFiResults = c_c'; 

    % % Load and store pertubations/ RVs. Size: n_samp x n_stochastic parameters
    load rv.mat
    xi_low = rv; 

    % I think this is appropriate. 

    xi_ref = rv; 
    % load and store high fidelity results (coarse and fine)
    load c_f.mat
    highFiResults = c_f'; 

    % % Grid points of interest - not certain on this one. 
    % load grid_coarse.mat % 1 x 50 grid points (x coordinates... ) 
    load grid_coarse.mat
    x_l = grid_c;

    % fine grid points 1 x 314 points
    load grid_fine.mat
    x_h = grid_f;

    gridpt_l = 1:length(x_l); 
    gridpt_h = 1:length(x_h); 

    % u_ref should be number of samples x spatial/temporal dimension of the
    % problem. size: n_samples x n_gridpoints
    u_ref = highFiResults(gridpt_h,:)';

    % u_low should be number of samples x spatial/temporal dimension of the
    % problem. size: n_samples x n_gridpoints
    u_low = lowFiResults(gridpt_l,:)';

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Gas Turbine
if run_GT == 1
  
    % load and store low fidelity results
    % two options for now, 1/21 cost (assembledRunMid_40_2) or 1/51 (assembledRunMid_40) cost. 
    
    lowFiResults = importdata('assembledRunMid_40_2'); 
    lowFiResults = lowFiResults';

    % % Load and store pertubations/ RVs. Size: n_samp x n_stochastic parameters
    y_samp = load('uniform_8.mat');
    xi_low = y_samp.uniform_data;

    
    y_samp = load('uniform_16_high.mat');
    xi_ref = y_samp.uniform_data;
    
    % load and store high fidelity results (coarse and fine)

    highFiResults = importdata('assembledRunMid_110_high'); 
    highFiResults = highFiResults';

%     load fine coordinates and reorder both data and coordinates by y. 
    x_h = importdata('assembledCoordsMid_110_high');
    [~,idx_h] = sort(x_h(2,:));
    x_h = x_h(2,idx_h); 
    highFiResults = highFiResults(idx_h,:); 
    
    x_l = importdata('assembledCoordsMid_40_2');
    [~,idx_l] = sort(x_l(2,:));
    x_l = x_l(2,idx_l);
    lowFiResults = lowFiResults(idx_l,:); 
    
    
    %     data_base = data_base(:,idx); 

    gridpt_l = 1:length(x_l); 
    gridpt_h = 1:length(x_h); 

    % u_ref should be number of samples x spatial/temporal dimension of the
    % problem. size: n_samples x n_gridpoints
    u_ref = highFiResults(gridpt_h,:)';

    % u_low should be number of samples x spatial/temporal dimension of the
    % problem. size: n_samples x n_gridpoints
    u_low = lowFiResults(gridpt_l,:)';
    

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desired polynomial order of PCE

p = 4; 

% 2 for LIB
% 5 for airfoil, 4 for LDC
% sets range of number of hi-fidelity samples. Can also just look at a
% single N_hi
% % N_hi = [10 30 50 70 90 110];
N_hi = [20]; 
% N_hi = [10, 15, 20, 25];

% Q: What exactly is the approximation rank? 
% A: this is the truncation rank of the approximation. Corresponds to the
% number of eigenvalues/eigenvectors retained

% set approximation rank of bi-fidelity model - can be a vector of values
% or a single rank. To see low performance of the bi-fidelity model set to
% something really low like r =2;

% r = [20];

% r = [1 3 8];
r = [5]; 
% r = [1 3 4 5 6]; 

% tolerance on residual used in spgl1 solver
sigma = .001;

% number of repetitions of psi_low the experiment for a given r and N_hi. When you
% want to run a full blown experiment set to 100. 
n_reps = 2; 

%% PCE stuff

% index_pc = nD_polynomial_array
% assembles index of polynomials for a given stochastic dimension and 
% polynomail order. 

% one column for each variable (ie dimension d). Number of rows for set of
% P basis functions based on order pol (or p) ie 0 0, 1 0, 0 1 for p = 1
% yielding P = 3

index_pc = nD_polynomial_array(size(xi_ref,2),p); 

% size of PC basis (set of P basis functions 
P = size(index_pc,1);

t_setup = toc; 
% disp('completed setup, time: ') 
fprintf('Setup with time : %d s.\n', t_setup);
%% Reference Solution

% measurement matrix used by the reference model
% size: number of data points x number of basis functions

psi_ref = zeros(size(xi_ref,1),P);

for isim=1:size(xi_ref,1)
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi_ref(isim,:),index_pc);
    psi_ref(isim,:) = crow_ref(1:P);
end


% Q: verbosity?
% A: detail of output messages

% Important parameters for spgl1 toolbox
% verbosity = 0 suppresses spgl1 command line output
opts = spgSetParms('iterations',10000,'verbosity',0);

% Reference PCE coefficents via l_1 minimization
% spg_mmv is a function in the spgl1 toolbox -Solve multi-measurement basis pursuit denoise

% perhaps change this to spgl1.. 

c_ref = spg_mmv(psi_ref,u_ref,sigma*norm(u_ref),opts);

% replace mmv with this and comment out transpose

% weights = get_matrix_weights(psi_ref);
% Psiweights = psi_ref*weights;
% 
% c_ref = zeros(size(u_ref,2), P);
% sigma_vec = zeros(size(u_ref,2),1); 
% 
% 
% % Step through each coordinate point
% for i_points=1:size(u_ref,2)
% % % sigma is truncation error of PCE (approximated)
% sigma_1 =  cross_val_sigma(psi_ref,u_ref(:,i_points));
% c_ref(i_points,:) = weights*spg_bpdn(Psiweights,u_ref(:,i_points),sigma_1*norm(u_ref(:,i_points)),opts);
% sigma_vec(i_points) = sigma_1; 
% end

% c_ref = psi_ref\u_ref; 
% transposed for convenience
c_ref = c_ref';

t_ref = toc - t_setup; 
fprintf('Reference solution : %d s.\n', t_ref);

%% Low-fidelity solution 

% measurement matrix used by low-fidelity model
psi_low = zeros(size(xi_low,1),P);

for isim=1:size(xi_low,1)
    % piset evaluates legendre polynomials at xi.
    crow_low = piset(xi_low(isim,:),index_pc);
    psi_low(isim,:) = crow_low(1:P);
end

% spgl1 options
opts = spgSetParms('iterations',10000,'verbosity',0);

% low fidelity coefficients solved via l_1 minimization
c_low = spg_mmv(psi_low, u_low, sigma*norm(u_low),opts);

% c_low = psi_low\u_low; 

% transposed for convenience
c_low = c_low';

% replace mmv with this and comment out transpose

% weights = get_matrix_weights(psi_low);
% Psiweights = psi_low*weights;
% 
% c_low = zeros(size(u_low,2), P);
% sigma_vec = zeros(size(xi_low,1),1); 
% % Step through each coordinate point
% 
% for i_points=1:size(u_low,2)
% % % sigma is truncation error of PCE (approximated)
% sigma_1 =  cross_val_sigma(psi_low,u_low(:,i_points));
% c_low(i_points,:) = weights*spg_bpdn(Psiweights,u_low(:,i_points),sigma_1*norm(u_low(:,i_points)),opts);
% sigma_vec(i_points) = sigma_1; 
% end


t_low = toc - t_ref; 
fprintf('Low fidelity solution : %d s.\n', t_low);
%% Relative errors of reference and low fidelity solution

% mean of reference solution
mean_ref = c_ref(:,1); 

% mean of low fidelity solution
mean_low = c_low(:,1);

% variance of ref solution
% Q: why 2:end?
% A: first term c_ref(1) is the mean. Including this would give the first
% moment. To adjust to central moment subtract mean^2 or leave out first
% term. 
var_ref = sum(c_ref(:,2:end).^2,2); 
std_ref = sqrt(sum(c_ref(:,2:end).^2,2)); 

% variance of low fidelity solution
var_low = sum(c_low(:,2:end).^2,2); 
std_low = sqrt(sum(c_low(:,2:end).^2,2)); 

% relative error of the mean as predicted by low fidelity model w.r.t.
% reference PCE


% relative error of the mean as predicted by low fidelity model w.r.t ref
% PCE. Identify grid points. 

% appears to be specific to cavity problem

%%%%%%%%%%%%%%%%%%%%%%% LDC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_LDC == 1

% Interpolate low error]
n_cell_l = 8;
x_l = linspace(0,1,n_cell_l+1);
x_l = 0.5*(cos(pi*(x_l-1)/2)+1);

n_cell_h = 32;
x_h = linspace(0,1,n_cell_h+1);
x_h = 0.5*(cos(pi*(x_h-1)/2)+1);

mean_low_int = interp1q(x_l', mean_low, x_h');

% mean_low_err = norm(mean_low - mean_ref(1:4:end))/norm(mean_ref); 
mean_low_err = norm(mean_low_int - mean_ref)/norm(mean_ref); 

% relative error of the variance as predicted by low fidelity model w.r.t
% ref PCE
% var_low_err = norm(var_low - var_ref(1:4:end))/norm(var_ref); 
var_low_int = interp1q(x_l', var_low, x_h');

var_low_err = norm(var_low_int - var_ref)/norm(var_ref); 
end

%%%%%%%%%%%%%%%%%%%%%%% LIB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_LIB == 1
    
    % interpolate 

    mean_low_int = interp1(x_l', mean_low, x_h', 'linear', 'extrap');
    var_low_int = interp1(x_l', var_low, x_h', 'linear', 'extrap');

    mean_low_err = norm(mean_low_int - mean_ref)/norm(mean_ref);
    % relative error of the variance as predicted by low fidelity model w.r.t.
    % reference PCE
    var_low_err = norm(var_low_int - var_ref)/norm(var_ref);

    t_errors = toc - t_low; 
    fprintf('relative errors : %d s.\n', t_errors);

end

%%%%%%%%%%%%%%%%%%%%%%% airfoil  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_airfoil ==1 
    
    mean_low_err = norm(mean_low - mean_ref)/norm(mean_ref); 
    var_low_err = norm(var_low - var_ref)/norm(var_ref); 
    
end

%%%%%%%%%%%%%%%%%%%%%%% gas turbine  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if run_GT ==1 
 
    mean_low_int = interp1(x_l', mean_low, x_h', 'linear', 'extrap');
    var_low_int = interp1(x_l', var_low, x_h', 'linear', 'extrap');

    mean_low_err = norm(mean_low_int - mean_ref)/norm(mean_ref);
    % relative error of the variance as predicted by low fidelity model w.r.t.
    % reference PCE
    var_low_err = norm(var_low_int - var_ref)/norm(var_ref);

    t_errors = toc - t_low; 
    fprintf('relative errors : %d s.\n', t_errors);
    
end

%% Bi-fidelity and High fidelity solutions

% consider storing as a vector and reshaping... 
% means as a function of r
mean_mean_bi_err = zeros(length(r)*length(N_hi),1); 
mean_mean_hi_err = zeros(length(r)*length(N_hi),1); 



% mean_mean_bi_err = zeros(length(r), length(N_hi)); 
% mean_mean_hi_err = zeros(length(r), length(N_hi)); 

var_mean_bi_err = zeros(length(r), length(N_hi)); 
var_mean_hi_err = zeros(length(r), length(N_hi)); 

mean_var_bi_err = zeros(length(r), length(N_hi)); 
mean_var_hi_err = zeros(length(r), length(N_hi)); 

var_var_bi_err = zeros(length(r), length(N_hi)); 
var_var_hi_err = zeros(length(r), length(N_hi)); 

mean_lam_hi = cell(length(r), length(N_hi)); 
mean_lam_ref = cell(length(r), length(N_hi)); 
mean_lam_low = cell(length(r), length(N_hi)); 
        
disp('about to start r and N loop') 

for i_rank = 1:length(r)        % index range of rank of bi-fidi_hi
    for i_hi = 1:length(N_hi)   % index range of number of high fidelity simulations
     
%         initialize parameters
%         bi-fidelity mean at grid points, generated ntotal times
        mean_bi_n1 = zeros(size(u_ref,2), n_reps); 
        var_bi_n1 = zeros(size(u_ref, 2), n_reps); 
        mean_hi_n1 = zeros(size(u_ref, 2), n_reps); 
        var_hi_n1 = zeros(size(u_ref, 2), n_reps); 

        mean_bi_n0 = zeros(size(u_low,2), n_reps); 
        var_bi_n0 = zeros(size(u_low, 2), n_reps); 
        mean_hi_n0 = zeros(size(u_low, 2), n_reps); 
        var_hi_n0 = zeros(size(u_low, 2), n_reps); 
        
        
%         eigenvalues of the reference, high and low fidelity models (also
%         low estimate)
        lam_ref1 = zeros(r(i_rank), n_reps); 
        lam_hi1 = zeros(r(i_rank), n_reps); 
        lam_low1 = zeros(r(i_rank), n_reps); 
        lam_ref0 = zeros(r(i_rank), n_reps); 
        lam_hi0 = zeros(r(i_rank), n_reps);        
        lam_low0 = zeros(r(i_rank), n_reps); 
        
%         bi-fid error in predicting the mean w.r.t ref mean
        mean_bi_err = zeros(1, n_reps); 
        
%         hi-fid error in predicting the mean w.r.t ref mean
        mean_hi_err = zeros(1, n_reps); 
        
%         bi-fid error in predicting variance w.r.t ref
        var_bi_err = zeros(1, n_reps); 
        
%         hi-fid error in predicting variance w.r.t ref
        var_hi_err = zeros(1, n_reps); 
      
        % Store something for error bound here... 
        
% %         sv_rp1_low_vec1 = zeros(ntotal,1); % singular value r+1 (or k+1)
% %         sv_rp0_low_vec1 = zeros(ntotal,1); % singular value r (or k)
        
        
%         index for number of times you are repeating the experiment.
%         ntotal is defined in the first section of this script
        
        % Bifidelity stuff:
%         H_mat = zeros(size(u_ref,2),N_hi);
        
        
        for i_rep = 1:n_reps    %index number of repetitions
            
%             randomly select sample set out of all available high fidelity
%             samples
% use randperm rather than datasample? 

            
            sample = datasample(1:size(u_ref,1), N_hi(i_hi), 'Replace', false); 
%             sample = 1:N_hi; 
            
        
            
            %% High fidelity model
%             random variables at sample
            xi_hi = xi_ref(sample,:); 
            
%             high fidelity model evaluations at sample
            u_hi = u_ref(sample,:); 
            
%             measurement matrix at sample - used to form high fid PCE
            psi_hi= psi_ref(sample,:); 
            
%             spgl1 optionsxi_low
            opts = spgSetParms('iterations',10000,'verbosity',0);
            
%             Hi fid PCE solved for via ell_1 minimization
            c_hi = spg_mmv(psi_hi,u_hi,sigma*norm(u_hi),opts);
            
%             c_hi = psi_hi\u_hi; 

%             transposed for convenience
            c_hi = c_hi';
                      
%             weights = get_matrix_weights(psi_hi);
%             Psiweights = psi_hi*weights;
% 
%             c_hi = zeros(size(u_hi,2), P);
%             sigma_vec = zeros(length(u_hi(1,:)),1); 
%             % Step through each coordinate point
%             for i_points=1:size(u_hi,2)
%             % % sigma is truncation error of PCE (approximated)
%             sigma_1 =  cross_val_sigma(psi_hi,u_hi(:,i_points));
%             c_hi(i_points,:) = weights*spg_bpdn(Psiweights,u_hi(:,i_points),sigma_1*norm(u_hi(:,i_points)),opts);
%             sigma_vec(i_points) = sigma_1; 
%             end

%             mean of high fidelity model
            mean_hi_n1(:,i_rep) = c_hi(:,1);
            
%             variance of high fidelity model
            var_hi_n1(:,i_rep) = (sum(c_hi(:,2:end).^2,2));
            std_high = sqrt((sum(c_hi(:,2:end).^2,2))); 
            
%             2nd moment (subtract mean^2 to find variance) - equivilent to
%             excluding first term
            
            %% Bi fidelity model
%             Outputs:
%             mean_bi_n1 - mean of the bi fidelity model at repitition n
%             var_bi_n1 - variance of the bi fid model at repitition n
%             c_bi1 - bi-fid model new reduced basis. Not saved here
%             lam_ref1 - eigenvalues of hi-fid model
%             lam_hi1 - eigenvalues of the hi-fid model
%             lam_low1 -eigenvalues of low-fid model
            
%        xi_high     Inputs:
%             u_hi - high-fid samples at sample
%             psi_hi - hi-fid measurement matrix
%             hi-fid - PCE coefficients
%             low-fid - PCE coefficients
%             reference PCE coefficients
%             current rank of approximation (index i_rank)
%             P - possibly no longer used
%             sigma - tolerance on l_2 residual
            
            % for purposes of error bound, generate bi-fidelity estimate of
            % low fidelity data. ]
            % Use all low fidelity samples. Error arrises from truncation.
% %              [mean_bi_n0(:,i_rep),var_bi_n0(:,i_rep),c_bi0,...
% %                 lam_ref0(:,i_rep),lam_hi0(:,i_rep),lam_low0(:,i_rep), psi_bi_low0,sv_rp1_low_vec1(i_rep),sv_rp0_low_vec0(i_rep)]= ...
% %                 BR_FN(u_low,psi_low,c_low,c_low,c_ref,P,r(i_rank),sigma);      
             [mean_bi_n0(:,i_rep),var_bi_n0(:,i_rep),c_bi0,...
                lam_ref0(:,i_rep),lam_hi0(:,i_rep),lam_low0(:,i_rep), psi_bi_low0]= ...
                BR_FN(u_low,psi_low,c_low,c_low,c_ref,P,r(i_rank),sigma); 
            
            % Now I have c_bi0. and psi_bi_low0, use to approximate L


            % generate bi-fidelty estimate based on based on many low-fidelity samples
%             and some high fidelity. This is the general approach
% %             [mean_bi_n1(:,i_rep),var_bi_n1(:,i_rep),c_bi1,...
% %                 lam_ref1(:,i_rep),lam_hi1(:,i_rep),lam_low1(:,i_rep), psi_bi_low,sv_rp1_low_vec1(i_rep),sv_rp0_low_vec1(i_rep),alpha2]= ... 
% %                 BR_FN(u_hi,psi_hi,c_hi,c_low,c_ref,P,r(i_rank),sigma);   
            [mean_bi_n1(:,i_rep),var_bi_n1(:,i_rep),c_bi1,...
                lam_ref1(:,i_rep),lam_hi1(:,i_rep),lam_low1(:,i_rep), psi_bi_low,alpha2]= ... 
                BR_FN(u_hi,psi_hi,c_hi,c_low,c_ref,P,r(i_rank),sigma);               
            % generate bi-fidelty based on selelct hi-fidelity samples. (c_low
%             changed to c_hi) % as would have access to in practice. 
            
            [mean_bi_n2,var_bi_n2,c_bi2,...
                lam_ref2,lam_hi2,lam_low2, psi_bi_hi,alpha2]= ...
                BR_FN(u_hi,psi_hi,c_hi,c_hi,c_ref,P,r(i_rank),sigma);

            % generate bi-fidelty based on ref samples - gold standard. 
            % would not have access in practice
            
            [mean_bi_n3,var_bi_n3,c_bi3,...
                lam_ref3,lam_hi3,lam_low3, psi_bi_ref,alpha2]= ...
                BR_FN(u_hi,psi_hi,c_hi,c_ref,c_ref,P,r(i_rank),sigma);

%           store errors of mean and variance calculations
            
            mean_bi_err(i_rep) = norm(mean_bi_n3 - mean_ref)/norm(mean_ref); 
            mean_hi_err(i_rep) = norm(mean_hi_n1(:,i_rep) - mean_ref)/norm(mean_ref); 
            var_bi_err(i_rep) = norm(var_bi_n1(:,i_rep)- var_ref)/norm(var_ref); 
            var_hi_err(i_rep) = norm(var_hi_n1(:, i_rep) - var_ref)/norm(var_ref); 
            
            std_bi_err = norm(sqrt(var_bi_n1(:,i_rep))-sqrt(var_ref))/norm(sqrt(var_ref)); 
            mean_bi_err;
     
            fprintf('repetition : %d of %d \n', i_rep, n_reps);
            fprintf('Timer : %d s.\n', toc);
            
 1;
% % % %             %% Error bound stuff
% % % %             % make script that does error estimate or some such. 
% % % %             % Store stuff for errors... 
% % % %             % Need to look at H and L and Lhat
% % % %                         
% % % %             % Calculate true error?
% % % %             H_mat_n = u_hi'; % size M x n (just 20 samples for instance)
% % % %             H_mat = u_ref';
% % % %             
% % % %             % Bi fidelity estimate of all high fidelity samples. 
% % % %             % First construct bi fidelity basis: 
% % % %             psi_bi_est = psi_ref(:,2:end)*alpha2';
% % % %             psi_bi_est = [ones(size(psi_ref, 1), 1) psi_bi_est]; 
% % % %             
% % % %             % remove r+1 column before solving for coefficients.
% % % %             psi_bi_est = psi_bi_est(:,1:end-1); 
% % % %                         
% % % %             H_mat_bi = c_bi1*psi_bi_est';
% % % %             
% % % %             1;
% % % %             
% % % %             L_mat = u_low'; % size m x N (ie all 200 samples)
% % % %             L_mat_n =  u_low(sample,:)'; 
% % % %     
% % % %             L_mat_bi = c_bi0*psi_bi_low0';
% % % %             
% % % % %             % What should I take svd of? 
% % % % % %             not psi... as done currently in BR_FN. 
% % % % %             % svd of L in practical error bounds paper. Try this for now. 
% % % % %             S = svd(L_mat); 
% % % % %             sv_rp0_low_vec1(i_rep) = S(r);
% % % % %             sv_rp1_low_vec1(i_rep) = S(r+1);
% % % %             
% % % %             1;
% % % %             
% % % %             % Calculate epsilon(tau)
% % % %             c = 1; % should this be N_low/N_hi or account for r? 
% % % % 
% % % %             exp1 = [-5:0.1:5]; 
% % % %             tau_vec = 10.^exp1';
% % % % 
% % % %             epsilon = zeros(length(tau_vec), 1); 
% % % %             bound_1_vec = zeros(length(tau_vec), 1);
% % % %             E_vec = zeros(length(tau_vec), 1); 
% % % %             T_vec = zeros(length(tau_vec), 1); 
% % % %             
% % % %             L = norm((L_mat-L_mat_bi),2);
% % % %             % use all pairs of epsilon tau to find minimum bound
% % % % 
% % % % % %             for i_tau = 1:length(tau_vec)
% % % % % %                 epsilon(i_tau) = c*max(eig((H_mat_n'*H_mat_n) -tau_vec(i_tau)*(L_mat_n'*L_mat_n))); 
% % % % % %                 
% % % % % %                 E_vec(i_tau) = sqrt(tau_vec(i_tau)*sv_rp1_low_vec1(i_rep)^2+epsilon(i_tau));
% % % % % %                 T_vec(i_tau) = sqrt(tau_vec(i_tau) + epsilon(i_tau)/sv_rp0_low_vec1(i_rep)^2);
% % % % % %                 
% % % % % %                 bound_1_vec(i_tau) = 2*E_vec(i_tau)+T_vec(i_tau)*L;
% % % % % %             end
% % % % % %             
% % % % % %             [bound_1,index_1] = min(bound_1_vec); 
% % % %             
% % % %             % bound_1 is aspects of bound related to truncation
% % % %             
% % % %             % id
% % % % %             figure 
% % % % %             loglog(tau_vec, epsilon)
% % % % % 
% % % % %             figure
% % % % %             hold on
% % % % %             plot(E_vec,'r')
% % % % %             plot(T_vec,'k')
% % % % %             hold off          
% % % % %          % Now caclulate parts accounted for by sampling error. 
% % % %             
% % % %             
% % % %             1;
% % % %             
% % % % %             Calculate true Error
% % % %         error_true = norm(H_mat- H_mat_bi, 2);
            
           
        end
        % Store errors as a function of rank r.
        
        1; 
        
        mean_mean_bi_err(i_rank,i_hi) = mean(mean_bi_err);
        mean_mean_hi_err(i_rank,i_hi) = mean(mean_hi_err);
        
        var_mean_bi_err(i_rank,i_hi) = var(mean_bi_err);
        var_mean_hi_err(i_rank,i_hi) = var(mean_hi_err);
        
        mean_var_bi_err(i_rank,i_hi) = mean(var_bi_err);
        mean_var_hi_err(i_rank,i_hi) = mean(var_hi_err);
        
        var_var_bi_err(i_rank,i_hi) = var(var_bi_err);
        var_var_hi_err(i_rank,i_hi) = var(var_hi_err);
    
%         average eigenvalues at a given rank. 
%         Note: saved as a cell because the number of eigenvalues changes 
%         w.r.t rank
        mean_lam_hi{i_rank,i_hi} = num2cell(mean(lam_hi1,2)); 
        mean_lam_ref{i_rank,i_hi} = num2cell(mean(lam_ref1,2)); 
        mean_lam_low{i_rank,i_hi} = num2cell(mean(lam_low1,2)); 
        
%         1; 
        fprintf('N_h  : %d \n', N_hi(i_hi));
        fprintf('r : %d s.\n', r(i_rank));
%         fprintf('Timer : %d s.\n', toc);
    end 
%         Save for a given approximation rank and number of hi-fid samples
%         N_hi
%         save(['BR_Nhi_' num2str(N_hi(i)) '_r_' num2str(r(k))]);
end
    
%     save for a given approximation rank
%     save(['BR_r_' num2str(r(k))])

% %     Print outputs 
%     mean_var_bi_err
%     mean_var_hi_err
%     var_low_err
%     mean_mean_bi_err
%     mean_mean_hi_err
%     mean_low_err

% ^ should fix this

% saves everything
save('Results/BR_FN','mean_mean_bi_err' , 'mean_mean_hi_err' ,...
    'var_mean_bi_err', 'var_mean_hi_err', ...
    'mean_var_bi_err' , 'mean_var_hi_err' ,...
    'var_var_bi_err', 'var_var_hi_err', ...
    'mean_lam_hi', 'mean_lam_ref', 'mean_lam_low','N_hi',...
    'var_low_err','mean_low_err','psi_bi_low', 'psi_bi_hi')

% Plot eigenvalues 
cov_low = c_low(:, 2:end) * c_low(:, 2:end)'; 
[v_low, lam_low1] = eigs(cov_low, 8); 
lam_low1 = abs(diag(lam_low1)); 

figure
semilogy(lam_low1/max(lam_low1), '-ok', 'Linewidth', 2);
set(gca,'Fontsize',14);
xlabel('index i', 'interpreter', 'latex', 'fontsize', 24);
ylabel('$\lambda_i$', 'interpreter', 'latex', 'fontsize', 24);

mean_mean_hi_err
mean_mean_bi_err


% % save results for time snapshot plot
mean_bi = mean(mean_bi_n3,2);
mean_hi = mean(mean_hi_n1,2);
var_bi = mean(var_bi_n1,2);
var_hi = mean(var_hi_n1,2);
% 
% mean_bi = mean_bi_n1(:,45);
% mean_hi = mean_hi_n1(:,45);
x_l;
x_h;

% if run_airfoil == 1
%     mean_low_int = mean_low;
%     var_low_int = var_low;
% end

if run_LDC == 1
    
end

mean_ref; 
mean_bi;
mean_hi;

figure
hold on
plot(x_h, mean_ref, 'k', 'linewidth', 2);
plot(x_l, mean_low, 'r', 'linewidth', 2);
plot(x_h, mean_bi, 'b', 'linewidth', 1);
plot(x_h, mean_hi, 'c', 'linewidth', 2);

figure
hold on
plot(x_h, var_ref, 'k', 'linewidth', 2);
plot(x_l, var_low, 'r', 'linewidth', 2);
plot(x_h, var_bi, 'b', 'linewidth', 1);
plot(x_h, var_hi, 'c', 'linewidth', 2);


% 
% figure
% hold on
% plot(x_l, u_low, 'r', 'linewidth', 2);
% plot(x_h, mean_hi_n1, 'b', 'linewidth', 2);
% plot(x_h, mean_bi_n1, 'c', 'linewidth', 1);

% u_low_samp = u_low(sample,:); 
% save('x_l','x_l')
% save('LIB_Nh_10_snapshot', 'x_h', 'mean_low_int', 'mean_ref', 'mean_bi','mean_hi','sample','u_low','u_hi','var_ref','var_low_int','var_bi','var_hi','u_ref', 'x_l')
% save('Airfoil_Nh_10_snapshot', 'x_h', 'mean_low_int', 'mean_ref', 'mean_bi','mean_hi','sample','u_low','u_hi')
% save('LDC_Nh_10_snapshot', 'x_h', 'mean_low_int', 'mean_ref', 'mean_bi','mean_hi','sample','u_low','u_hi')