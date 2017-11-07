function [E,E_ED,E_real]=Energy_X_RBM2(a,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)
%% Initialization
E=0;
RE=0;

E_ED=0;
RE_ED=0;

E_real=0;
RE_real=0;

N_par=N_up+N_dn;

%%
     for i2=1:2^(2*N_sites)
     x2(2*N_sites:-1:1,1)=bitget(i2-1,2*N_sites:-1:1);
     x2=(x2-0.5).*2;
     for j2=1:2^N_y %binary order
     y2(N_y:-1:1,1)=bitget(j2-1,N_y:-1:1);
     y2=(y2-0.5).*2;
         eN_up=zeros(N_sites,N_sites);
         eN_dn=zeros(N_sites,N_sites);
         for i=1:N_sites
             eN_up(i,i)=exp(a(i,1)*x2(i,1));
             eN_dn(i,i)=exp(-1*a(i+N_sites,1)*x2(i,1));
         end    
         
%          phi_1(:,1:N_up)=Phi_T(:,1:N_up).*exp(x1'*w*y1);
%          phi_1(:,N_up+1:N_par)=Phi_T(:,N_up+1:N_par).*exp(x1'*w*y1);
%          phi_1=Phi_T;
%          phi_2(:,1:N_up)=Phi_T(:,1:N_up).*exp(x2'*w*y2);
%          phi_2(:,N_up+1:N_par)=Phi_T(:,N_up+1:N_par).*exp(x2'*w*y2);

         phi_1=Phi_T;
         phi_2(:,1:N_up)=exp(x2'*w*y2)*eN_up*Proj_k*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par)=exp(x2'*w*y2)*eN_dn*Proj_k*phi_1(:,N_up+1:N_par);
         phi_2=Proj_k*phi_2;
         
        %% ED
         Phi_2=kron(phi_2(:,1:N_up),phi_2(:,N_up+1:N_par));
         Phi_1=kron(phi_1(:,1:N_up),phi_1(:,N_up+1:N_par));
         
         % H_ED
         H_K=zeros(N_sites,N_sites);
         H=zeros(2*N_sites,2*N_sites);
         H_K(1,2)=-2;
         H_K(2,1)=-2;
         H=kron(H_K,H_K);
         H(1,1)=U;
         H(2*N_sites,2*N_sites)=U;
         
         [psi_nonint,E_nonint_m] = eig(H);
         E_nonint_v=diag(E_nonint_m);
         e_ed=E_nonint_v(1);
         re_ed=psi_nonint(:,1)'*psi_nonint(:,1);
         
         E_ED=E_ED+e_ed;
         RE_ED=RE_ED+re_ed;
         
         % real
         Phi_2_re=Phi_2/sqrt(Phi_2'*Phi_2);
         e_real=Phi_1'*H*Phi_2;
         re_real=Phi_1'*Phi_2;
         
         E_real=E_real+e_real;
         RE_real=RE_real+re_real;
         
        %%  calculate the single-particle Green's function matrix for each spin:
         re_up=det(phi_1(:,1:N_up)'*phi_2(:,1:N_up));
         re_dn=det(phi_1(:,N_up+1:N_par)'*phi_2(:,N_up+1:N_par));
         re=re_up*re_dn;
         
%          invO_matrix_up=inv(phi_1(:,1:N_up)'*phi_2(:,1:N_up));
%          invO_matrix_dn=inv(phi_1(:,N_up+1:N_par)'*phi_2(:,N_up+1:N_par));
         
         G_up=phi_2(:,1:N_up)*phi_1(:,1:N_up)';
         G_dn=phi_2(:,N_up+1:N_par)*phi_1(:,N_up+1:N_par)';

        %% calculate the potential energy:
         n_int=(diag(G_up)).'*diag(G_dn);
         potentialEnergy=n_int*U;

        %% calculate the kinetic energy:
         kineticEnergy=sum(sum(H_k.'.*(G_up+G_dn))); % note the element-wise multiplication

        %% calculate the total energy:
         e(i2,j2)=potentialEnergy+kineticEnergy;
         E=E+e(i2,j2);
         RE=RE+re;
         
     end
     end
 
 E=E/RE;
 E_ED=E_ED/RE_ED;
 E_real=E_real/RE_real;
 
 
end