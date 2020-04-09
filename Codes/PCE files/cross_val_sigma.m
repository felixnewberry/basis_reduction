function [sigma] = cross_val_sigma(A,b)
% This function attempts to find an optimal value of sigma for the bpdn
% problem: 
%    minimize ||x||_1 subject to ||A*x-b||_2 <= sigma
% Inputs:
% A: is an M x P sensing matrix, where M is the number of
% data samples in time and P is the number of basis functions.
% y: is an M x 1 data bector corrupted with noise.
% Output:
% sigma is an approximate of the optimal tolerance for the spg_bpdn problem

R = 4; %Number of folds or resampling
[M,P] = size(A); %Number of Data snapshots
N_sigmas = 20; %Number of sigma tolerances to attempt
opts = spgSetParms('iterations',2*P,'verbosity',0);
sigmas = linspace(0,5,N_sigmas);
sigmas = 10.^(-sigmas);
%sigmas = linspace(1e-10,1,N_sigmas); % Relative error Tolerances to attempt
%opts = spgSetParms('verbosity',0);  % Turn off the SPGL1 log output
errs = zeros(N_sigmas,1);

N_r = floor(0.8*M);
N_v = M - N_r;
for i = 1:N_sigmas
    sigma = sigmas(i);
    val_errs = zeros(R,1);
    for r = 1:R
%         Uncomment these lines, comment the next too to change from random
%         to fold validation. 

%         val_inds = (1+(r-1)*floor(R^(-1)*M)):r*floor(R^(-1)*M);
%         if r == R
%             val_inds = (1+(r-1)*floor(R^(-1)*M)):M;
%         end
%         recon_inds = setdiff(1:M,val_inds);
        recon_inds = datasample(1:M,N_r,'Replace',false);
        val_inds = setdiff(1:M,recon_inds);
        A_r = A(recon_inds,:);
        A_v = A(val_inds,:);
        b_r = b(recon_inds,:);
        b_v = b(val_inds,:);
        weights = get_matrix_weights(A_r);
        Xi = weights*spg_bpdn(A_r*weights,b_r,sigma*norm(b_r),opts);
        val_errs(r) = norm(A_v*Xi - b_v)/norm(b_v);
    end
    errs(i) = mean(val_errs);
end
% figure
% plot(sigmas,errs,'-o')

[~,min_ind] = min(errs);
sigma = errs(min_ind);

end

