% Generate LDC high-fidelity samples

% Save inputs_vec with delta_nu, delta_u and nx

% Generate an additional 1800 samples.. 



% delta_nu = 0; delta_u = 0; nx = 64; 
% save('LDC_data/inputs_vec_high', 'delta_nu', 'delta_u', 'nx')
% 
% % Save u_nu_vec_2.mat (name it high_1800)
% % nu_vec and u_top_vec
% 
% % nu_vec is 0.01 +- 10% 
% % u_top_vec is 1 +- 20% 
% 
% xi_hi = rand(2000,2); 
% nu_vec = (1+0.1*(2*xi_hi(:,1)-1))*0.01;
% u_top_vec = (1+0.2*(2*xi_hi(:,2)-1))*1;
% 
% save('LDC_data/u_nu_vec_hi', 'xi_hi', 'nu_vec', 'u_top_vec')

tic
% Call python/fenics
pyFi_ensemble = 'sudo python3.6 run_LDC_ensemble.py';
system(pyFi_ensemble); 
toc