function [err_E,err_T, err_1, err_2] = br_E_T_dev(B,A, Ir, psi_bi)
%Compute E and T to test whether this is much tighter

% Inputs: 
% B - low-fidelity samples
% A - highfidelity samples
% r - optimal rank that minimizes error - found in iteration over epsilon
% tau
% for comparing to using n samples I should maybe use n... 

% Outputs: 
% err_E
% err_T

%%% Following 'Practical Error Bounds', from lemma 3. in appendix. 

[~, ~, V] = svd(B); 

% find the T that minimizes the total error
T = A*V(:,1:Ir)*V(:,1:Ir)'*pinv(B); 

err_T = norm(T); 
err_E = norm(A-T*B); 

err_1 = norm(T*B*pinv(psi_bi)*psi_bi - A*pinv(psi_bi)*psi_bi);
err_2 = norm(T*B-A)*norm(pinv(psi_bi)*psi_bi); 
% norm(pinv(psi_bi)*psi_bi) should be 1


end

