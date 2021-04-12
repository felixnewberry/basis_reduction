function [c_hat] = SP( Psi, K, u , weight,cross_val_type)
% Subspace Persuit
% Inputs: Psi  NxP measurement matrix
%       : K = approximate bound on signal sparsity such that K >= s.
%         set K = 0 if bound is unknown to cross validate K.
%       : u is an Nx1 vector of measurements (QOI samples) 
%       : weight == true for normalized support estimates
%       : cross_val_type, set this to 'random' for random resampling, or
%         set this to 'fold' for k-fold cross-validation
% Output: c Px1 vector of PCE coefficients
%INITIALIZATION:   
c_hat = zeros(size(Psi,2),1);
stop = 0; num_iter=0; max_iter = 20*size(Psi,2);
if nargin < 5 
   cross_val_type = 'fold'; 
end
if nargin < 4 || weight == false
   W = eye(size(Psi,2)); 
   weight = false;
else 
   W = get_matrix_weights(Psi); 
   Psi =  Psi*W; 
end
if K == 0
    if strcmp('fold',cross_val_type)
        [K,~] = cross_val_SP_fold(Psi,u,weight);
    end
    if strcmp('random',cross_val_type)
        [K,~] = cross_val_SP(Psi,u,weight);
    end
end

%Initial support estimation
x = zeros(size(Psi,2),1);
corr = abs(Psi'*u); [vals,~] = sort(corr,'descend');
I = find(corr >= vals(K));
%Initial residual calculation
x(I) =  pinv(Psi(:,I))*u;  
ur0 = u - Psi*x;

while stop == 0 
    corr = abs(Psi'*ur0);  [vals,~] = sort(corr,'descend'); %corr = abs(Psi'*ur0);
    II = find(corr >= vals(K));II = union(I,II);
    %LSP iteration. 
    x = zeros(size(Psi,2),1);  x(II) = pinv(Psi(:,II))*u;  
    %Updated support estimation
    I0 = I;  [vals,~] = sort(abs(x),'descend');
    I =  find( abs(x) >= vals(K));  
    %Update Residual  
    x = zeros(size(Psi,2),1); x(I) = pinv(Psi(:,I))*u;
    ur = u - Psi*x; 
    num_iter=num_iter+1;
    %Check stopping criteria
    if norm(ur,2) >= norm(ur0,2)                      
        I = I0;                                          
        stop = 1;
    end
    if  num_iter == max_iter              
        stop = 1; 
    end
    ur0 = ur;
end
c_hat(I) = pinv(Psi(:,I))*u; %Alternatively Psi(:,I)\u;
c_hat = W*c_hat;

end