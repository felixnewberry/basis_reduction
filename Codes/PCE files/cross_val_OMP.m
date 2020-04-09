function [ TOL ] = cross_val_OMP( Psi, u, weight )
%Cross-validation function via OMP
% Inputs:
% Psi , NxP measurement matrix 
% u , Nx1 QoI vector

R = 4; %Number of 'folds'
N = size(Psi,1);
N_r = floor(0.8*N); %Numer of reconstruction samples
N_sigmas = 20; %Number of tolerances to attempt
sigmas = linspace(0,5,N_sigmas);
sigmas = 10.^(-sigmas);

OMP_errs = zeros(N_sigmas,1);
if nargin < 3 
    weight = false;
end

for i = 1:N_sigmas
    TOL = sigmas(i);
    val_errs = zeros(R,1);
    for r = 1:R
        recon_inds = datasample(1:N,N_r,'Replace',false);
        val_inds = setdiff(1:N,recon_inds);
%         if r == R
%             val_inds = (1+(r-1)*floor(R^(-1)*N)):N;
%         else
%             val_inds = (1+(r-1)*floor(R^(-1)*N)):r*floor(R^(-1)*N);
%         end
%         recon_inds = setdiff(1:N,val_inds);
        Psi_r = Psi(recon_inds,:);
        Psi_v = Psi(val_inds,:);
        u_r = u(recon_inds);
        u_v = u(val_inds);
        c_hat = OMP(Psi_r,TOL,u_r,weight);
        val_errs(r) = norm(Psi_v*c_hat - u_v)/norm(u_v);
    end
    OMP_errs(i) = mean(val_errs);
end
%figure
%loglog(sigmas,OMP_errs,'o-') 


[~,min_ind] = min(OMP_errs);
TOL =  OMP_errs(min_ind);

end

