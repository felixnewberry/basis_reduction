function [c, psi, err_val] = my_pce(xi, p, u, sigma, pc_solver, pc_val)
% Construct Polynomial Chaos Expansion
% Inputs: 
% xi        random variables
% p         polynomial order
% d         stochastic dimension
% u         sampled data - n_samps x n_points 
% sigma     spg solver tolerance
% pc_solver       use least squares, mmv or spg to solve (0, 1 and 2)
% pc_val    report validation error 1: 3/4 samples to train, 1/4 val. 0: no
% validation

% Outputs: 
% c         PCE coefficients
% psi       PCE basis
% err_val   Validation error

1; 

% if pc_val == 1
%     N_total = n_sim; 
%     N_test = round(0.75*N_total); 
%     N_val = N_total - N_test; 
% else
%     N_test = N_sim; 
% end

n_points = size(u,2); 
n_samps = size(xi,1); 

if pc_val == 1
    n_train = round(0.75*n_samps);
    n_val = n_samps-n_train; 
else
    n_train = n_samps; 
end

% Assemble index of polynomials
index_pc = nD_polynomial_array(size(xi,2),p); 


% size of PC basis (set of P basis functions 
P = size(index_pc,1);

% Construct polynomial basis
psi_train = zeros(n_train,P);

for isim=1:n_train
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi(isim,:),index_pc);
    psi_train(isim,:) = crow_ref(1:P);
end

opts = spgSetParms('iterations',10000,'verbosity',0);

if pc_solver == 0
    c = psi_train\u(1:n_train,:); 
    c = c'; 
elseif pc_solver == 1
    % Solve PCE coefficents via l_1 minimization
    c = spg_mmv(psi_train,u(1:n_train,:),sigma*norm(u(1:n_train,:)),opts);
    c = c';
else
    
    1;
    
    
    weights = get_matrix_weights(psi_train);
    Psiweights = psi_train*weights;

    
%     sigma_vec = zeros(n_points,1); 
    
    % Use parfor loop if stepping through individual spatial points
    % Can I use u and Psiweights in the parfor loop?
    norm_u_vec = vecnorm(u(1:n_train,:),2); 
    % Step through each coordinate point
    
    c_data{n_points} = []; 
    
    parfor i_points=1:n_points
        % % sigma is truncation error of PCE (approximated)
        
        opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);
        
%         sigma_1 =  cross_val_sigma(Psiweights,u(:,i_points));
        delta = sigma*norm_u_vec(i_points);
        c_data{i_points}.c_vec = weights*spg_bpdn(Psiweights,u(1:n_train,i_points),delta,opts);

%         sigma_vec(i_points) = sigma_1; 
%         sigma_1 =  cross_val_sigma(psi,u(:,i_points));
%         c(i_points,:) = weights*spg_bpdn(Psiweights,u(:,i_points),sigma_1*norm(u(:,i_points)),opts);
%         sigma_vec(i_points) = sigma_1; 
    end
    
    c = zeros(n_points, P);
    for i_points=1:n_points
        c(i_points,:) = c_data{i_points}.c_vec; 
    end
end

if pc_val == 1
    psi_val = zeros(n_val,P);
    for isim=1:n_val
        crow_val = piset(xi(isim+n_train,:),index_pc);
        psi_val(isim,:) = crow_val(1:P);
    end
    err_val = norm(u(n_train+1:end,:)-psi_val*c')/norm(u(n_train+1:end,:));
    psi = [psi_train; psi_val]; 
else
    psi = psi_train;
    err_val = NaN; 
end

end

