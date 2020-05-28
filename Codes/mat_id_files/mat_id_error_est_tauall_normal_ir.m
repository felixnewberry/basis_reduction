% function [delta_eps, ahat_error_est] = mat_id_error_est(B_R, A_R, normC, err_Bhat, sb)
% delta_eps = set of epsilon and corresponding delta for isometry conditions
% ahat_error_est = estimate of error
% rnk = rank of approximation to use
% B_R = Sub set of Low Fidelity Vectors
% A_R = Associated High fidelity Matrix
% normC = norm of coefficient matrix
% err_Bhat = error in bhat approximation
% sb = singular values of B

%Code updated from v7. Follow theory:
%  %  eps(t) = lam_max (A'A- delta*B'B)
%  % returns min_de1 = min_de2, such that term 1 and 2 are using same
%  points (t,eps)
% ir = returns index k of rho_k

function [delta_eps, ahat_error_bound, ir] = mat_id_error_est_tauall_normal_ir(B_R, A_R, normC, err_Bhat, sb,N,R)
    %Changed epsilon delta computation to a stopping criteria and inflation factor
    delta_min = .1; % 1 is a fine minimum usually
    delta_inflation = 0.01; % times larger each iteration
    max_delta_checks = 5*10^8; % Maximum number of delta to test Will also stop if epsilon increases

    % Generate epsilon-delta curve. This should be a considerably more stable computation.
    % I moved this outside of an inner loop where it was running once for each rank.
    % Step 1: delta = 0
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
    while count <= max_delta_checks && delta_eps(end,1)<10^8%&& delta_eps(end,2) < delta_eps(end-1,2)% stop if too many delta or epsilon becomes bad
        count = count + 1; % increment delta count
        delta_inc = delta_inflation*delta; % My best attempt to insure stability here
        delta = delta + delta_inc; % New Delta
        
        
        G2 = G - delta*BRMat; % New G
        G2 = (G2+G2')/2; % Enforce symmetry
        %[~, Lam] = function_lanczos_cov_rand(G2,size(B_R,2));
        Lam = eig(G2);
        Lam = diag(Lam);
        Lam = max(0,Lam);
        epsilon = max(diag(Lam));
        delta_eps = [delta_eps; delta epsilon];
    end

    delta_eps(:,2) = delta_eps(:,2).*N/R;
    num_eps_delta = size(delta_eps,1); % The number of (delta,epsilon) points computed
    
    
    r2 = size(sb,1)-1;
    r2 = min(size(sb,1)-1, rank(diag(sb)));
%     if r2>20
%         r2=20;
%     end

% % % % %UPDATE: for each point (tau, eps), optimize bound over r
    % For each r, minimize each term of bound contribution of (delta, epsilon(delta))
    ahat_error_bound = zeros(num_eps_delta,1);
    ir = zeros(num_eps_delta,1);
    bound = zeros(r2,2); 
    i_r = 1:(r2);
    i_rp = 2:(r2+1);
    for i_tau = 1:num_eps_delta
        tau = delta_eps(i_tau,1);
        ep =  delta_eps(i_tau,2);
        bound(:,1) = (1+normC)*sqrt((tau)*(sb(i_rp)).^2  + ep) ;
        bound(:,2) = sqrt( tau + ep./(sb(i_r).^2)) * err_Bhat;    
        [val , ind ] = min(abs(bound(:,1))+abs(bound(:,2))); % True minimizer 
        ahat_error_bound(i_tau) = val;
        ir(i_tau) = ind;
        
    end
    
    
end



