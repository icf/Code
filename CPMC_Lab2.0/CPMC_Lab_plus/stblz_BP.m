function Phi = stblz_BP(Phi, N_up, N_par)
% function [Phi, O] = stblz(Phi, N_wlk, O, N_up, N_par)
% Perform the modified Gram-Schmidt orthogonalization to stabilize the walkers
% Inputs:
%   Phi: the whole ensemble of walkers
%   N_wlk: the number of walkkers
%   O: the array of overlaps of all walkers
%   N_up: the number of spin up electrons
%   N_par: the total number of electrons
% Outputs:
%   Phi: the stabilized ensemble of walkers
%   O: the updated array of overlaps
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ?014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% Perform the QR decomposition on each walker
% Keep only the Q matrices and discard the R matrices
    [Phi(:,1:N_up),R_up]=qr(Phi(:,1:N_up),0);
    % for the spin down sector:
    [Phi(:,N_up+1:N_par),R_dn]=qr(Phi(:,N_up+1:N_par),0);
end