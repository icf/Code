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
     for y2=1:N_y %binary order
     
     n_phi=n_phi+1;
         eN_up=zeros(N_sites,N_sites);
         eN_dn=zeros(N_sites,N_sites);
         
         temp=a(y2,:);
         temp_T=temp';
         for i=1:N_sites
             eN_up(i,i)=exp(temp_T(i,1));
             eN_dn(i,i)=exp(-1*temp_T(i,1));
         end  
         
%          [psi_nonint,E_nonint_m] = eig(H_k);
%          Phi_T1=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn)); %Trivial K State.
%          phi_1(:,:,1)=Phi_T1;
%          
%          Phi_T2=horzcat(psi_nonint(:,1:N_up-1),psi_nonint(:,N_up+1),psi_nonint(:,1:N_dn)); %Trivial K State.
%          phi_1(:,:,2)=Phi_T2;
%       
%          Phi_T3=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn-1),psi_nonint(:,N_dn+1)); %Trivial K State.
%          phi_1(:,:,3)=Phi_T3;
%    
%          Phi_T4=horzcat(psi_nonint(:,1:N_up-1),psi_nonint(:,N_up+1),psi_nonint(:,1:N_dn-1),psi_nonint(:,N_dn+1)); %Trivial K State.
%          phi_1(:,:,4)=Phi_T4;
%          
%          Phi_T5=horzcat(psi_nonint(:,1:N_up-2),psi_nonint(:,N_up+1:N_up+2),psi_nonint(:,1:N_dn)); %Trivial K State.
%          phi_1(:,:,5)=Phi_T5;
%          
%          Phi_T6=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn-2),psi_nonint(:,N_dn+1:N_dn+2)); %Trivial K State.
%          phi_1(:,:,6)=Phi_T6;
%       
%          Phi_T7=horzcat(psi_nonint(:,1:N_up-2),psi_nonint(:,N_up+1:N_up+2),psi_nonint(:,1:N_dn-1),psi_nonint(:,N_dn+1)); %Trivial K State.
%          phi_1(:,:,7)=Phi_T7;
%    
%          Phi_T8=horzcat(psi_nonint(:,1:N_up-1),psi_nonint(:,N_up+1),psi_nonint(:,1:N_dn-2),psi_nonint(:,N_dn+1:N_dn+2)); %Trivial K State.
%          phi_1(:,:,8)=Phi_T8;
         
         phi_1=Phi_T;
         phi_2(:,1:N_up,n_phi)=w(y2,1)*eN_up(:,:)*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par,n_phi)=w(y2,1)*eN_dn(:,:)*phi_1(:,N_up+1:N_par);
     end

    %%  calculate the single-particle Green's function matrix for each spin:
    for i=1:n_phi
        for j=1:n_phi
            re_up=det(phi_2(:,1:N_up,i)'*phi_2(:,1:N_up,j));
            re_dn=det(phi_2(:,N_up+1:N_par,i)'*phi_2(:,N_up+1:N_par,j));
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
 E_real=E_real/RE_real;
 E_ED=E_ED/RE_ED;
 
end