function  Phi = HF_H_pin(n,H_k,N_sites,N_up,N_dn,U,Lx,Ly,Lz,kx,ky,kz,tx,ty,tz)
% Calculate Density at each site, icf 2017/10/16
%% Initialization
U2_diag=0;
U1_up=zeros(N_sites,N_sites);
U1_dn=zeros(N_sites,N_sites);
H_HF_up=zeros(N_sites,N_sites);
H_HF_dn=zeros(N_sites,N_sites);
psi_nonint_up=zeros(N_sites,N_sites);
psi_nonint_dn=zeros(N_sites,N_sites);

for i=1:N_sites
    U2_diag=U2_diag-0.5*U*n(i)*n(i+N_sites);
end

for i=1:N_sites
    U1_up(i,i)=U*n(i+N_sites); %+U2_diag/(N_up+N_dn);
    U1_dn(i,i)=U*n(i); %+U2_diag/(N_up+N_dn);
end

H_pin_up=H_PIN_up(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz);
H_pin_dn=H_PIN_up(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz);

H_HF_up=H_k+U1_up+H_pin_up;
H_HF_dn=H_k+U1_dn+H_pin_dn;

%% Calculation
[psi_nonint_up,D_up] = eig(H_HF_up);
[psi_nonint_dn,D_dn] = eig(H_HF_dn);
E_nonint_up=diag(D_up);
E_nonint_dn=diag(D_dn);
E=sum(E_nonint_up(1:N_up))+sum(E_nonint_dn(1:N_dn));


Phi=horzcat(psi_nonint_up(:,1:N_up),psi_nonint_dn(:,1:N_dn));

end