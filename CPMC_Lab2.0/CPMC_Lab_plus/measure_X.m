function [e, Proj] = measure_X(H_k, phi, Phi_T, Phi_T2, invO_matrix_up, invO_matrix_dn, invO_matrix_up2, invO_matrix_dn2, N_up, N_par, U)
% function e = measure(H_k, phi, Phi_T, invO_matrix_up, invO_matrix_dn, N_up, N_par, U)
% Calculate the mixed estimator for the ground state energy of a walker
% Inputs:
%   H_k: the one-body kinetic Hamiltonian
%   phi: the matrix of a single walker
%   Phi_T: the matrix of the trial wave function
%   invO_matrix_up: the inverse of the spin up sector of the walker's overlap matrix 
%   invO_matrix_dn: the inverse of the spin down sector of the walker's overlap matrix 
%   N_up: the number of spin up electrons
%   N_par: the total number of electrons of both spins
%   U: the on-site repulsion strength in the Hubbard model
% Outputs:
%   e: the mixed estimator for the ground state energy of the input walker
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ?014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

    %%  calculate the single-particle Green's function matrix for each spin:
    invO_matrix_up=inv(Phi_T(:,1:N_up)'*phi(:,1:N_up));
    invO_matrix_dn=inv(Phi_T(:,N_up+1:N_par)'*phi(:,N_up+1:N_par));
    temp_up=phi(:,1:N_up)*invO_matrix_up;
    temp_dn=phi(:,N_up+1:N_par)*invO_matrix_dn;
    G_up=temp_up*Phi_T(:,1:N_up)';
    G_dn=temp_dn*Phi_T(:,N_up+1:N_par)';
    %
    G_up2=Phi_T2(:,1:N_up)'*phi(:,1:N_up);
    G_dn2=Phi_T2(:,N_up+1:N_par)'*phi(:,N_up+1:N_par);
    %% Re
    re_up_T2=det(Phi_T2(:,1:N_up)'*phi(:,1:N_up));
    re_dn_T2=det(Phi_T2(:,N_up+1:N_par)'*phi(:,N_up+1:N_par));
    re_T2=re_up_T2*re_dn_T2;
    %
    re_up_T=det(Phi_T(:,1:N_up)'*phi(:,1:N_up));
    re_dn_T=det(Phi_T(:,N_up+1:N_par)'*phi(:,N_up+1:N_par));
    re_T=re_up_T*re_dn_T;
    %%
    Proj=det(G_up2)*det(G_dn2)/re_T;

    %% calculate the potential energy:
    n_int=(diag(G_up)).'*diag(G_dn);
    potentialEnergy=n_int*U;

    %% calculate the kinetic energy:
    kineticEnergy=sum(sum(H_k.'.*(G_up+G_dn))); % note the element-wise multiplication

    %% calculate the total energy:
    e=potentialEnergy+kineticEnergy;

end
