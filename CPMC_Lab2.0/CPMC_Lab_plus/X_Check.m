function [E,N]=X_Check(Lx,Ly,Lz,kx, ky,kz, tx, ty,tz,N_up,N_dn,U,N_y,w_trace,phi_2)
%% Initialization
H_k=H_K(Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);

N_sites=Lx*Ly*Lz;
N_par=N_up+N_dn;

w(:,1)=w_trace(end,:,1);
n_phi=N_y;

E=0;
N=0;
RE=0;

for i=1:n_phi
    for j=1:n_phi
            re_up=w(i,1)*w(j,1)*det(phi_2(:,1:N_up,i)'*phi_2(:,1:N_up,j));
            re_dn=w(i,1)*w(j,1)*det(phi_2(:,N_up+1:N_par,i)'*phi_2(:,N_up+1:N_par,j));
            re=re_up*re_dn;
            
            invO_matrix_up=inv(phi_2(:,1:N_up,i)'*phi_2(:,1:N_up,j));
            invO_matrix_dn=inv(phi_2(:,N_up+1:N_par,i)'*phi_2(:,N_up+1:N_par,j));
         
           %%  calculate the single-particle Green's function matrix for each spin:
            temp_up=phi_2(:,1:N_up,j)*invO_matrix_up;
            temp_dn=phi_2(:,N_up+1:N_par,j)*invO_matrix_dn;
         
            G_up=temp_up*phi_2(:,1:N_up,i)';
            G_dn=temp_dn*phi_2(:,N_up+1:N_par,i)';
            G_up=G_up*re;
            G_dn=G_dn*re;

          %% calculate the potential energy:
           n_int=(diag(G_up)).'*diag(G_dn);
           potentialEnergy=n_int*U/re;

          %% calculate the spin density
           n_up=abs(diag(G_up));
           n_dn=abs(diag(G_dn));
           n(1:N_sites)=n_up(1:N_sites);
           n(N_sites+1:2*N_sites)=n_dn(1:N_sites);
 
          %% calculate the kinetic energy:
           kineticEnergy=sum(sum(H_k.'.*(G_up+G_dn))); % note the element-wise multiplication

          %% calculate the total energy:
           e(i,j)=potentialEnergy+kineticEnergy;
           
           N=N+n;
           E=E+e(i,j);
           RE=RE+re;
         
        end
    end
 
 N=N/RE;
 E=E/RE;
     
end