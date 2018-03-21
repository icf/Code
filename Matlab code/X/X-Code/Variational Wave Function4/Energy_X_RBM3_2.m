function [E,E_ED,E_real,N]=Energy_X_RBM3_2(a,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)
%% Initialization
E=0;
RE=0;

E_ED=0;
RE_ED=0;

E_real=0;
RE_real=0;

N=0;

n_phi=0;

N_par=N_up+N_dn;

%%
     for j2=1:N_y %binary order
%          if j2<=2^(N_y-1)
%              y2=zeros(2^(N_y-1),1);
%              y2(j2,1)=1;
%          else
%              y2=zeros(2^(N_y-1),1);
%              y2(j2-2^(N_y-1),1)=1;
%              if j2-2^(N_y-1)+1 > 2^(N_y-1)
%                 y2(1,1)=1;
%              else
%                 y2(j2-2^(N_y-1)+1,1)=1;
%              end
%          end
         y2(N_y:-1:1,1)=bitget(j2-1,N_y:-1:1);
         y2=(y2-0.5);
             
         n_phi=n_phi+1;
         eN_up=zeros(N_sites,N_sites);
         eN_dn=zeros(N_sites,N_sites);
         
         temp=y2'*a;
         temp_T=temp';
         for i=1:N_sites
             eN_up(i,i)=exp(temp_T(i,1));
             eN_dn(i,i)=exp(-1*temp_T(i,1));
         end  
         
         phi_1=Phi_T;
         phi_2(:,1:N_up,n_phi)=eN_up(:,:)*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par,n_phi)=eN_dn(:,:)*phi_1(:,N_up+1:N_par);
     end

    %%  calculate the single-particle Green's function matrix for each spin:
    for i=1:n_phi
%             re_up=det(phi_2(:,1:N_up,i)'*phi_2(:,1:N_up,i));
%             re_dn=det(phi_2(:,N_up+1:N_par,i)'*phi_2(:,N_up+1:N_par,i));
%             re=re_up*re_dn;
            
            invO_matrix_up=inv(phi_2(:,1:N_up,i)'*phi_2(:,1:N_up,i));
            invO_matrix_dn=inv(phi_2(:,N_up+1:N_par,i)'*phi_2(:,N_up+1:N_par,i));
         
           %%  calculate the single-particle Green's function matrix for each spin:
            temp_up=phi_2(:,1:N_up,i)*invO_matrix_up;
            temp_dn=phi_2(:,N_up+1:N_par,i)*invO_matrix_dn;
         
            G_up=temp_up*phi_2(:,1:N_up,i)';
            G_dn=temp_dn*phi_2(:,N_up+1:N_par,i)';
            G_up=G_up;
            G_dn=G_dn;

          %% calculate the potential energy:
           n_int=(diag(G_up)).'*diag(G_dn);
           potentialEnergy=n_int*U;

          %% calculate the spin density
           n_up=abs(diag(G_up))*abs(w(i,1));
           n_dn=abs(diag(G_dn)*abs(w(i,1)));
           n(1:N_sites)=n_up(1:N_sites);
           n(N_sites+1:2*N_sites)=n_dn(1:N_sites);
 
          %% calculate the kinetic energy:
           kineticEnergy=sum(sum(H_k.'.*(G_up+G_dn))); % note the element-wise multiplication

          %% calculate the total energy:
           e(i,1)=(potentialEnergy+kineticEnergy)*abs(w(i,1));
           
           N=N+n;
           E=E+e(i,1);
           RE=RE+1*abs(w(i,1));
         
    end
 
 N=N/RE;
 E=E/RE;
 E_real=E_real/RE_real;
 E_ED=E_ED/RE_ED;
 
end


