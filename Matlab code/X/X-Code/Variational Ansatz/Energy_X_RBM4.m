function [E,E_ED,E_real]=Energy_X_RBM4(a,w,b,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)
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
     for i2=1:2^(N_sites)
     x2(N_sites:-1:1,1)=bitget(i2-1,N_sites:-1:1);
     x2=(x2-0.5).*2;
     for j2=1:2^N_y %binary order
     y2(N_y:-1:1,1)=bitget(j2-1,N_y:-1:1);
     y2=(y2-0.5).*2;
     n_phi=n_phi+1;
         eN_up=zeros(N_sites,N_sites);
         eN_dn=zeros(N_sites,N_sites);
         for i=1:N_sites
%              temp_up=0;
%              temp_dn=0;
%              for j=1:N_sites
%                  temp_up=temp_up+a(i,j)*x2(j,1);
%                  temp_dn=temp_dn-a(i,j)*x2(j,1);
%              end
             eN_up(i,i)=exp(a(i,1)*x2(i,1)+b(i,1));
             eN_dn(i,i)=exp(-a(i,1)*x2(i,1)+b(i+N_sites,1));
         end    

         phi_1=Phi_T;
         phi_2(:,1:N_up,n_phi)=exp(x2'*w*y2)*eN_up*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par,n_phi)=exp(x2'*w*y2)*eN_dn*phi_1(:,N_up+1:N_par);
         phi_2(:,:,n_phi)=phi_2(:,:,n_phi);
         
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