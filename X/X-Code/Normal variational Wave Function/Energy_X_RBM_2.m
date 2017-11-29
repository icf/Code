function [E,E_ED,E_real]=Energy_X_RBM_2(a,b,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)
%% Initialization
E=0;
RE=0;

E_ED=0;
RE_ED=0;

E_real=0;
RE_real=0;

n_phi=0;

N_par=N_up+N_dn;

%%
     for j2=1:2^N_y %binary order
     y2(N_y:-1:1,1)=bitget(j2-1,N_y:-1:1);
     y2=(y2-0.5).*2;
     
     n_phi=n_phi+1;
         eN_up=zeros(N_sites,N_sites,N_up);
         eN_dn=zeros(N_sites,N_sites,N_dn);
         
         temp=y2'*a;
         temp_T=temp';
         for i=1:N_sites
             for j=1:N_up
                 eN_up(i,i,j)=exp(temp_T(i+(j-1)*N_sites,1));
             end
             for j=1:N_dn
                 eN_dn(i,i,j)=exp(-1*temp_T(i+(j-1)*N_sites,1));
             end
         end    
         
         phi_1=Phi_T;
         for j=1:N_up
             phi_2(:,j,n_phi)=exp(y2'*w)*eN_up(:,:,j)*phi_1(:,j);
         end
         for j=1:N_dn
             phi_2(:,N_up+j,n_phi)=exp(y2'*w)*eN_dn(:,:,j)*phi_1(:,N_up+j);
         end
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
 
          %% calculate the kinetic energy:
           kineticEnergy=sum(sum(H_k.'.*(G_up+G_dn))); % note the element-wise multiplication

          %% calculate the total energy:
           e(i,j)=potentialEnergy+kineticEnergy;
           
           E=E+e(i,j);
           RE=RE+re;
           
         %% ED
%           Phi_j=kron(phi_2(:,1:N_up,j),phi_2(:,N_up+1:N_par,j));
%           Phi_i=kron(phi_2(:,1:N_up,i),phi_2(:,N_up+1:N_par,i));
%          
%           % H_ED
%           H_K=zeros(N_sites,N_sites);
%           I=zeros(N_sites,N_sites);
%           H=zeros(2*N_sites,2*N_sites);
%           V=zeros(2*N_sites,2*N_sites);
%           I(1,1)=1;
%           I(2,2)=1;
%           H_K(1,2)=-2;
%           H_K(2,1)=-2;
%           H=kron(H_K,I)+kron(I,H_K);
%           H(1,1)=U;
%           H(2*N_sites,2*N_sites)=U;
%           V(1,1)=U;
%           V(2*N_sites,2*N_sites)=U;
% %           % CC_ED
% %           CC11=zeros(2*N_sites,2*N_sites);
% %           C_T=zeros(2*N_sites,2*N_sites);
% %           C_T(1,1)=1;
% %           C_T(2,2)=1;
% %           CC11=Phi_i'*C_T*Phi_j
%          
%          [psi_nonint,E_nonint_m] = eig(H);
%          E_nonint_v=diag(E_nonint_m);
%          e_ed=E_nonint_v(1);
%          re_ed=psi_nonint(:,1)'*psi_nonint(:,1);
%          
%          E_ED=E_ED+e_ed;
%          RE_ED=RE_ED+re_ed;
%          
%          % real
%          e_real=Phi_i'*H*Phi_j;
%          re_real=Phi_i'*Phi_j;
%          
%          E_real=E_real+e_real;
%          RE_real=RE_real+re_real;
         
        end
    end
 
 E=E/RE;
 E_real=E_real/RE_real;
 E_ED=E_ED/RE_ED;
 
end