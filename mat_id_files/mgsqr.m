function [Q,R,piv,k] = mgsqr(A,tol)
% Pivoted QR factorization via modified Gram-Schmidt.  Adapted from Golub
% and Van Loan, "Matrix Computations".
%
% Input: 
%  A   = low-rank matrix of size m-by-n
%  tol = stop when the largest column norm is below this value
%
% Output:
%  Q   = m-by-k matrix with orthonormal columns, s.t., Q^T Q = I
%  R   = k-by-n upper triangular matrix
%  piv = indices of pivot columns
%  k   = rank of the decomposition

  [m,n] = size(A);
  maxrnk = min(m,n);
  Q = zeros(m,maxrnk);
  R = zeros(maxrnk,n);
  
  piv = 1:n;
  c = sum(A.^2,1);
  [cmax,imax] = max(c);
  cmax = cmax(1);  imax = imax(1);

  for k = 1:maxrnk
    
    %fprintf('k = %d, imax = %d, cmax = %g\n', k, imax, cmax)
    
    % perform pivot
    tmp = piv(k);  piv(k) = piv(imax);  piv(imax) = tmp;
    tmp = A(:,k);  A(:,k) = A(:,imax);  A(:,imax) = tmp;
    tmp = R(:,k);  R(:,k) = R(:,imax);  R(:,imax) = tmp;
    tmp = c(k);    c(k) = c(imax);      c(imax) = tmp;
    
    % gram-schmidt sweep
    R(k,k) = norm(A(:,k));
    Q(:,k) = A(:,k) / R(k,k);
    R(k,(k+1):n) = Q(:,k)' * A(:,(k+1):n);
    for j = (k+1):n
      A(:,j) = A(:,j) - R(k,j) * Q(:,k);
    end
    
    % find next pivot 
    if k < maxrnk
      c((k+1):n) = sum(A(:,(k+1):n).^2,1);
      [cmax,imax] = max(c((k+1):n));
      cmax = cmax(1);  imax = imax(1);  imax = k + imax(1);
    end
    
    % check for convergence
    if abs(cmax) < tol
      Q = Q(:,1:k);
      R = R(1:k,:);
      return
    end
    
    
  end

end