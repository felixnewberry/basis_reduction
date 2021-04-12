function [A] = mat_normalize(A)

for i = 1:size(A,2)
    A(:,i) = A(:,i)./norm(A(:,i));
end

% D = diag(1./sqrt(diag(A)));
% A = D*A*D;

end

