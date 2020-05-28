function [P,ix] = matrixID(Y,tol)
% Compute the skeletonization of a matrix Y, s.t., Y = Y(:,ix)*P
% See demo.m for usage.

  % compute the QR factorization of Y
  [~,R,ix,k] = mgsqr(Y,tol);
  Pr = eye(size(Y,2));
  Pr = Pr(:,ix);
  ix = ix(1:k);

  % Non-lazy way to solve for T
  [Ru,Rs,Rv] = svd(R(1:k,1:k));
  rnk = min(size(Rs));
  Rs = Rs(1:rnk,1:rnk);
  Rs = Rs + Rs(1,1) * eps * eye(rnk);
  Rpinv = Rv(:,1:rnk) * diag(1./diag(Rs)) * Ru(:,1:rnk)';
  T = Rpinv * R(1:k,(k+1):end);
  
  % Lazy way to solve for T
  %T = R(1:k,1:k) \ R(1:k,(k+1):end);

  P = [eye(k),T];
  P = P * Pr';  
end