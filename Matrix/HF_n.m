function  n = HF_n(Phi,N_sites,N_up,N_dn,N_par)
% Calculate Density at each site, icf 2017/10/16
%% Initialization
invO_matrix_up=inv(Phi(:,1:N_up)'*Phi(:,1:N_up));
invO_matrix_dn=inv(Phi(:,N_up+1:N_par)'*Phi(:,N_up+1:N_par));
temp_up=Phi(:,1:N_up)*invO_matrix_up;
temp_dn=Phi(:,N_up+1:N_par)*invO_matrix_dn;
G_up=temp_up*Phi(:,1:N_up)';
G_dn=temp_dn*Phi(:,N_up+1:N_par)';

%% Calculation
n_up=abs(diag(G_up));
n_dn=abs(diag(G_dn));
n(1:N_sites)=n_up(1:N_sites);
n(N_sites+1:2*N_sites)=n_dn(1:N_sites);

end