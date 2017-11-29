function phi = stepwlk_BP(phi, x, N_sites, Proj_k_half, N_up, N_par, aux_fld)
% Perform one step of the random walk
% Inputs:
%   phi: the whole ensemble of walkers
%   N_wlk: the number of walkers
%   N_sites: the total number of lattice sites
%   w: the array of weights of all the walkers
%   O: the array of overlaps of all the walkers
%   E: the total energy of all walkers
%   W: the total weight of all walkers
%   H_k: the one-body kinetic Hamiltonian
%   Proj_k_half: the matrix of the operator exp(-deltau*K/2)
%   flag_mea: the flag (1 or 0) that specifies whether the energy should the measured in this step
%   Phi_T: the matrix of the trial wave function
%   N_up: the number of spin up electrons
%   N_par: the total number of electrons
%   U: the on-site repulsion strength in the Hubbard model
%   fac_norm: the exponent of the pre-factor exp(-deltau*(H-E_T))
%   aux_fld: the 2x2 matrix containing all the possible values of the quantity exp(gamma*s(sigma)*x_i) (used in V.m only)
% Outputs:
%   phi: the ensemble of walkers after propagation
%   w: the new array of weights of all walkers
%   O: the new array of overlaps of all walkers
%   E: the new total energy of all walkers
%   W: the new total weight of all walkers
%   
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ?014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% Propagate each walker:
    % propagate by the kinetic term exp(-1/2*deltau*K)
    phi=Proj_k_half*phi;
    % propagate each lattice site of a walker by the potential term:
    for j_site=1:N_sites
         x_spin=x(j_site);
         phi(j_site,1:N_up)=phi(j_site,1:N_up)*aux_fld(1,x_spin);
         phi(j_site,N_up+1:N_par)=phi(j_site,N_up+1:N_par)*aux_fld(2,x_spin); 
    end
    phi=Proj_k_half*phi;
            
end