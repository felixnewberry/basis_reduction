function [err_bound, ahat_error_est, err_bound_N, err_bound_n, err_T, ...
    err_E, omega_term, q_N, q_n, Ir] = ...
    br_bound(B_R, A_R, A_n, err_Bhat, sb, R, N, n, psi_bi)
%Compute the error bound for the basis reduction algorithm

% Inputs: 
% B_R = Sub set of Low Fidelity Vectors (R cols of mxN matrix B)
% A_R = Associated High fidelity Matrix  (R cols of MxN matrix A)
% A_n = High fidelity Matrix for omega term  (n cols of MxN matrix A)
% *** note, in practice we normalize matrices s.t. 
% A = U_hf/norm(U_hf,fro) and B = U_lf/norm(U_lf,fro), where U_hf and U_lf are the LF and HF data matrices
% err_Bhat = error in bhat approximation
% sb = singular values of B
% R = subset of data from which bound eps(tau) is estimated (R << N)
% N = total number of data columns
% n - number of samples used in bi-fidelity estimate (n << N)
% psi_bi = bi-fidelity basis of n samples to find coefficients

% Outputs: 
% e_bound = Estimated rror bound
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
[ahat_error_est,Ir] = min(ahat_error_bound);
% Components
err_E = sqrt((delta_eps(Ide(Ir),1))*(sb(Ir+1))^2  + delta_eps(Ide(Ir),2)); 
err_T = sqrt( delta_eps(Ide(Ir),1) + delta_eps(Ide(Ir),2)./(sb(Ir))^2); 

1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LS approximation with n < N samples 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n - number of samples used in estimate
% N - Total number of sims to estimate
% omega - take max for each point as surrogate
% q - solve from coherence - bit of a handful. 

% q: compute from condition 1.4 in "On the stability of Least Squares" Cohen et al. 
% Condition __ in our paper. 
% Desire the largest q (r in cohen) that satisfies the condition. 

% first - coherence, then solve via... 
% Should be a function to do for both n and N

% I need to change this psi_bi - do when I put it in a function.
% psi_bi_2 = 
% mu = max(sum(abs(psi_bi).^2,2)); 
% is this right? - should be orthonormal basis

1; 

% My method:
% A1 = mat_normalize(psi_bi(:,2:end)');
A1 = psi_bi(:,2:end)'; 
mu = max(sum(abs(A1).^2,1)); 

% Don't want to normalize. Test this. 

% This code does more, not sure that it's useful:
% [mu,mu_A_av] = coherence(psi_bi');

% Think I may be computing the coherence wrong? It looks like it should be
% this simple... 

% q_n is negative. Thoughts: 
% Hmm... 

% % does log in paper mean log10 or log e? 
% log e according to correction. Also matches intuition 

% Natural log
q_N = 0.5*(N*(3*log(3/2)-1)/(mu*log(N))-2); 
q_n = 0.5*(n*(3*log(3/2)-1)/(mu*log(n))-2);


% q_n = 0.5*(n*(3*log(3/2)-1)/(log(n))-2); 

1; 

% % Should the coherence be normalized?
% % psi_bi - L in paper. Arbitrary basis for Vm. For their analysis they
% % suppose that the basis is orthonormal. 
% % each basis vector is unit and orthogonal with the others. 
% A1 = mat_normalize(psi_bi');
% A2 = psi_bi'; 
% 
% mu1 = max(sum(abs(A1).^2,1)); 
% mu2 = max(sum(abs(A2).^2,1)); 
% 
% q_n1 = 0.5*(n*(3*log(3/2)-1)/(mu1*log(n))-2) %; 
% q_n2 = 0.5*(n*(3*log(3/2)-1)/(mu2*log(n))-2) %; 
% 
% % Should I maybe not be including the first column? - the ones
% % Test this: 
% A1 = mat_normalize(psi_bi(:,2:end)');
% A2 = psi_bi(:,2:end)'; 
% 
% mu1 = max(sum(abs(A1).^2,1)); 
% mu2 = max(sum(abs(A2).^2,1)); 
% 
% % test normal: 
% figure 
% hist(sum(abs(A1).^2,1)); 
% hold on
% hist(sum(abs(A2).^2,1)); 
% 
% 
% q_n1 = 0.5*(n*(3*log(3/2)-1)/(mu1*log(n))-2) %; 
% q_n2 = 0.5*(n*(3*log(3/2)-1)/(mu2*log(n))-2) %; 
% % mu2 is a result of a pretty significant outlier

% % should psi_bi be what was used with the n samples... consider. 
% n/log(n)
% 
% I should put a check here on sign
if q_N <= 0
    fprintf("Error: For N = %d, q_N < 0, (q_n = %6.2f), and condition is invalid. \n", N, q_N);
end
if q_n <= 0
    fprintf("Error: For n = %d, q_N < 0, (q_n = %6.2f), and condition is invalid. \n", n, q_n);
end

% Display the n and the q. 

1; 

% omega - maximum of each sample - estimate with H_n as surrogate
omega_term = sum(max(A_n,[],2).^4)^0.5;

%%% all together

err_bound = ahat_error_est+ 8*N^0.5*omega_term*(N^(-q_N)+n^(-q_n));
err_bound_N = 8*N^0.5*omega_term*(N^(-q_N)); 
err_bound_n = 8*N^0.5*omega_term*(n^(-q_n)); 
end

