% function [delta_eps, ahat_error_est,Ir, min_de1,min_de2] = mat_id_error_est(B_R, A_R, normC, err_Bhat, sb)
% delta_eps = [tau, eps] =set of epsilon and corresponding delta for isometry conditions
% ahat_error_est = estimate of error
% Ir - index of rank k (not the bi-fid rank) which minimized bound
% min_de1 = min_de2 = (tau, eps) that minimizes bound. Originally 2 points
% rnk = rank of approximation to use
% B_R = Sub set of Low Fidelity Vectors (R cols of mxN matrix B)
% A_R = Associated High fidelity Matrix  (R cols of MxN matrix A)
% *** note, in practice we normalize matrices s.t. 
% A = U_hf/norm(U_hf,fro) and B = U_lf/norm(U_lf,fro), where U_hf and U_lf are the LF and HF data matrices
% with all N columns. 
% normC = norm of coefficient matrix C_lf s.t. 
% U_lf \approx U_lf^{column skel}*
% C_lf
% err_Bhat = error in bhat approximation
% sb = singular values of B
% N = total number of data columns
% R = subset of data (n in the paper) from which bound eps(tau) is
% estimated

%Code updated from v7. Follow theory:
%  %  eps(t) = lam_max (A'A- delta*B'B)
%  % returns min_de1 = min_de2, such that term 1 and 2 are using same
%  points (t,eps)

function [delta_eps, ahat_error_est,Ir, min_de1,min_de2] = mat_id_error_est_one_normal(B_R, A_R, normC, err_Bhat, sb,N,R)
    %Changed epsilon delta computation to a stopping criteria and inflation factor
    delta_min = .1; % 1 is a fine minimum usually
    delta_inflation = 0.01; % times larger each iteration
    max_delta_checks = 5*10^5; % Maximum number of delta to test Will also stop if epsilon increases
%     max_delta_checks = 5*10^8; % Maximum number of delta to test Will also stop if epsilon increases

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
%     if r2>20
%         r2=20;
%     end

    % For each r, minimize each term of bound contribution of (delta, epsilon(delta))
    ahat_error_bound = zeros(r2,1);
    bound = zeros(num_eps_delta,2); 
    
    m_de1 = zeros(r2,2);
    m_de2 = zeros(r2,2);
    
    for i_r = 1:r2
        bound(:,1) = (1+normC)*sqrt((delta_eps(:,1))*(sb(i_r+1))^2  + delta_eps(:,2)) ;
        bound(:,2) = sqrt( delta_eps(:,1) + delta_eps(:,2)./(sb(i_r))^2) * err_Bhat;    
        ahat_error_bound(i_r) = min(abs(bound(:,1))+abs(bound(:,2))); % True minimizer
        %ahat_error_bound(i_r) = min(abs(bound(:,1)))+min(abs(bound(:,2))); % True minimizer
        
        [~,Ide1] =  min((abs(bound(:,1))+abs(bound(:,2))));
        [~,Ide2] =  min((abs(bound(:,1))+abs(bound(:,2))));
        
        m_de1(i_r) = Ide1;
        m_de2(i_r) = Ide2;
    end
    ahat_error_est = min(ahat_error_bound);
    [~,Ir] = min(ahat_error_bound);
    
    min_de1 = delta_eps(m_de1(Ir),:);
    min_de2 = delta_eps(m_de2(Ir),:);
    %get min de1 and de2
end



