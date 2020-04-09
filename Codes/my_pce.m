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


% Assemble index of polynomials
index_pc = nD_polynomial_array(size(xi,2),p); 

% size of PC basis (set of P basis functions 
P = size(index_pc,1);

% Construct polynomial basis
psi = zeros(size(xi,1),P);

for isim=1:size(xi,1)
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
    weights = get_matrix_weights(psi);
    Psiweights = psi*weights;

    c = zeros(size(u,2), P);
    sigma_vec = zeros(size(u,2),1); 
    % Step through each coordinate point
    for i_points=1:size(u,2)
        % % sigma is truncation error of PCE (approximated)
        
        opts = spgSetParms('iterations',8000,'verbosity',0,'optTol',1e-9,'bpTol',1e-9);
        delta = sigma*norm(u_samp(1:N,1));
        
%         sigma_1 =  cross_val_sigma(psi,u(:,i_points));
        c(i_points,:) = weights*spg_bpdn(Psiweights,u(:,i_points),delta,opts);
        sigma_vec(i_points) = sigma_1; 
%         sigma_1 =  cross_val_sigma(psi,u(:,i_points));
%         c(i_points,:) = weights*spg_bpdn(Psiweights,u(:,i_points),sigma_1*norm(u(:,i_points)),opts);
%         sigma_vec(i_points) = sigma_1; 
    end

end
end

