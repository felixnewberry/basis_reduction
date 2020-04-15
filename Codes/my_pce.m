function [c,psi] = my_pce(xi, p, u, sigma, pc_solver)
% Construct Polynomial Chaos Expansion
% Inputs: 
% xi        random variables
% p         polynomial order
% d         stochastic dimension
% u         sampled data
% sigma     spg solver tolerance
% pc_solver       use least squares, mmv or spg to solve (0, 1 and 2)

% Outputs: 
% c         PCE coefficients
% psi       PCE basis


n_points = size(u,2); 
n_samps = size(xi,1); 

% Assemble index of polynomials
index_pc = nD_polynomial_array(size(xi,2),p); 

% size of PC basis (set of P basis functions 
P = size(index_pc,1);

% Construct polynomial basis
psi = zeros(n_samps,P);

for isim=1:n_samps
%     piset evaluates a multi dimensional pc basis at xi. (legendre 
%     appropriate for  uniform RV expansion)
    crow_ref = piset(xi(isim,:),index_pc);
    psi(isim,:) = crow_ref(1:P);
end

opts = spgSetParms('iterations',10000,'verbosity',0);

if pc_solver == 0
    c = psi\u; 
elseif pc_solver == 1
    % Solve PCE coefficents via l_1 minimization
    c = spg_mmv(psi,u,sigma*norm(u),opts);
    c = c';
else
    
    1;
    
    
    weights = get_matrix_weights(psi);
    Psiweights = psi*weights;

    
%     sigma_vec = zeros(n_points,1); 
    
    % Use parfor loop if stepping through individual spatial points
    % Can I use u and Psiweights in the parfor loop?
    norm_u_vec = vecnorm(u,2); 
    % Step through each coordinate point
    
    c_data{n_points} = []; 
    
    parfor i_points=1:n_points
        % % sigma is truncation error of PCE (approximated)
        
        opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);
        
%         sigma_1 =  cross_val_sigma(Psiweights,u(:,i_points));
        delta = sigma*norm_u_vec(i_points);
        c_data{i_points}.c_vec = weights*spg_bpdn(Psiweights,u(:,i_points),delta,opts);

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
end

