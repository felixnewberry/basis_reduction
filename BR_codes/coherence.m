function [mu_A,mu_A_av] = coherence(A)
% 
mu = [];
A = mat_normalize(A);
G = A'*A;
%G = mat_normalize(G);

M = size(A,1);
L = size(A,2);
mu_bar = sqrt((L-M)/(M*(L-1)));

the_sum = 0;
num_inds = 0;

for i = 1:size(G,1)
    for j = 1:size(G,1)
        
       if i ~= j
          mu = [mu ; abs(A(:,i)'*A(:,j))/(norm(A(:,i))*norm(A(:,j)))];
       end
       
       if  mu_bar <= abs(G(i,j))
          the_sum = the_sum + abs(G(i,j));
          num_inds = num_inds+1;
       end
    end
        
end

mu_A = max(mu);

mu_A_av = the_sum/num_inds;

1;

    
end

