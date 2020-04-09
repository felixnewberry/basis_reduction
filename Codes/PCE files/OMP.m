function [ c_hat ] = OMP(Psi,TOL,u,weight)
% Standard Orthogonal Matching Persuit (OMP)
% Inputs: Psi  NxP measurement matrix
%       : u    Nx1 measured signal
%       : tol real number used to determine stopping criteria
%       : max_iter is the maximum number of iterations (use s)

if nargin < 4 || weight == false
   W = eye(size(Psi,2)); 
   weight = false;
else 
   W = get_matrix_weights(Psi); 
   Psi =  Psi*W; 
end

if TOL == -inf
    TOL = cross_val_OMP(Psi,u,weight);
end
P = size(Psi,2);
c_hat = zeros(P,1);
x = c_hat;
I = [];
stop = 0; iter = 1; max_iter = 2*size(Psi,2); % was 20 not 2.
ur = u;
while stop == 0
    corr = Psi'*ur;
    [~,i] =  max(abs(corr));
    I = union(I,i); %Merge support
    x = zeros(P,1); x(I) = pinv(Psi(:,I))*u;
    ur = u - Psi*x;
    iter = iter + 1;
    if TOL >= norm(ur,2)/norm(u) || iter == max_iter                                       
       stop = 1;
    end
end
c_hat(I) = pinv(Psi(:,I))*u;
c_hat = W*c_hat;
end

