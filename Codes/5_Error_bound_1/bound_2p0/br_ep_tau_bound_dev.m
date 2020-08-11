function [err_ep_tau_bound, err_E, err_T, Ir] = ...
    br_ep_tau_bound_dev(B_R, A_R, err_Bhat, sb, R, N)
%Compute the error bound for the basis reduction algorithm

% Inputs: 
% B_R = Sub set of Low Fidelity Vectors (R cols of mxN matrix B)
% A_R = Associated High fidelity Matrix  (R cols of MxN matrix A)

% *** note, in practice we normalize matrices s.t. 
% A = U_hf/norm(U_hf,fro) and B = U_lf/norm(U_lf,fro), where U_hf and U_lf are the LF and HF data matrices
% err_Bhat = error in bhat approximation
% sb = singular values of B
% R = subset of data from which bound eps(tau) is estimated (R << N)
% N = total number of data columns

% Outputs: 
% e_bound = Estimated error bound
% ahat_error_est = component assuming all N samples avaialable
% omega_term = max values of data
% q_n = term to compute error from n << N samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Epsilon Tau terms - assumes N samples available 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_min = .1; % 1 is a fine minimum usually
delta_inflation = 0.01; % times larger each iteration
max_delta_checks = 5*10^5; % Maximum number of delta to test Will also stop if epsilon increases

G = A_R'*A_R ;
G = (G+G')/2; % Enforce symmetry
G2 = G; % We reuse G, using G2 for intermediate computations

Lam = eig(G2);
Lam = diag(Lam); % non-absolute valued lambdas don't help
BRMat = B_R'*B_R; % Compares with Lam
BRMat = (BRMat+BRMat')/2; % Enfore Symmetry
[V_B, Lam_B] = eig(BRMat);
Lam_B = diag(Lam_B);
Lam_B = max(0,real(Lam_B)); % Enforce positive semi definiteness
BRMat = V_B*diag(Lam_B)*V_B';

epsilon = max(diag(Lam));
delta_eps = [0 epsilon]; % epsilon for delta = 0.

delta = delta_min;
G2 = G - delta*BRMat;
G2 = (G2+G2')/2; % Enforce symmetry
Lam= eig(G2);
Lam = max(0,diag(Lam)); % Remove all negative bits
epsilon = max(diag(Lam)); % Quick evaluation of the norm when we do it this way
delta_eps = [delta_eps; delta_min epsilon]; %#ok<*AGROW> % epsilon for minimum delta

count = 2; % Two values computed so far
while count <= max_delta_checks && delta_eps(end,2) < delta_eps(end-1,2)% stop if too many delta or epsilon becomes bad
    count = count + 1; % increment delta count
    delta_inc = delta_inflation*delta; % My best attempt to insure stability here
    delta = delta + delta_inc; % New Delta

    G2 = G - delta*BRMat; % New G
    G2 = (G2+G2')/2; % Enforce symmetry
    Lam = eig(G2);
    Lam = diag(Lam);
    Lam = max(0,Lam);
    epsilon = max(diag(Lam));
    delta_eps = [delta_eps; delta epsilon];
end

delta_eps(:,2) = delta_eps(:,2).*N/R;
num_eps_delta = size(delta_eps,1); % The number of (delta,epsilon) points computed

r2 = size(sb,1)-1;

% For each r, minimize each term of bound contribution of (delta, epsilon(delta))
ahat_error_bound = zeros(r2,1);
bound = zeros(num_eps_delta,2); 

Ide = zeros(r2,1);
% m_de2 = zeros(r2,2);

for i_r = 1:r2
%     bound(:,1) = (1+normC)*sqrt((delta_eps(:,1))*(sb(i_r+1))^2  + delta_eps(:,2)) ;
    bound(:,1) = 2*sqrt((delta_eps(:,1))*(sb(i_r+1))^2  + delta_eps(:,2)) ;

%     bound(:,2) = sqrt( delta_eps(:,1) + delta_eps(:,2)./(sb(i_r))^2) * err_Bhat;    
    bound(:,2) = sqrt( delta_eps(:,1) + delta_eps(:,2)./(sb(i_r))^2) * err_Bhat;    

    [ahat_error_bound(i_r), Ide(i_r)] = min(abs(bound(:,1))+abs(bound(:,2))); % True minimizer

end
% Total
[err_ep_tau_bound,Ir] = min(ahat_error_bound);

% Components
err_E = sqrt((delta_eps(Ide(Ir),1))*(sb(Ir+1))^2  + delta_eps(Ide(Ir),2)); 
err_T = sqrt( delta_eps(Ide(Ir),1) + delta_eps(Ide(Ir),2)./(sb(Ir))^2); 

end

