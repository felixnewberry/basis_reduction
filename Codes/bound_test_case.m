% Test Error Bound

% Problem setup

clear all
close all
clc

% What do I need access to? 

% Epsilon tau               Approximate with Grammians of L and H. Check
% this method. 

% sigma_{k+1}, sigma_k      Compute from L

% L                         Readily available

% \hat{L}                   Compute using bi-fidelity interpolate that is
% constructed with L_N - all the low-fidelity samples. 

% Theta                     Truncation error of KL - compute from
% eigenvalues of covariance matrix. (either during truncation or at some
% other point. Would be better to have entirely seperate. 

% mu -                      Coherence - copmute from bi-fidelity basis. 

%%% Inputs: 

% L, H_n, r, p (tune p as a preliminary step?)
% See how that goes. 