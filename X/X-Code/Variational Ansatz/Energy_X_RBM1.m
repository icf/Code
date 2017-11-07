function E=Energy_X_RBM1(a,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)
%% Initialization
E=0;
RE=0;
N_par=N_up+N_dn;

%%
     for i2=1:2^(2*N_sites)
     x2(2*N_sites:-1:1,1)=bitget(i2-1,2*N_sites:-1:1);
     x2=(x2-0.5).*2;
     for j2=1:2^N_y %binary order
     y2(N_y:-1:1,1)=bitget(j2-1,N_y:-1:1);
     y2=(y2-0.5).*2;
         eN2=zeros(2*N_sites,2*N_sites);
         for i=1:2*N_sites
             eN2(i,i)=exp(-1*a(i,1)*x2(i,1));
         end    
         
%          phi_1(:,1:N_up)=Phi_T(:,1:N_up).*exp(x1'*w*y1);
%          phi_1(:,N_up+1:N_par)=Phi_T(:,N_up+1:N_par).*exp(x1'*w*y1);
%          phi_1=Phi_T;
%          phi_2(:,1:N_up)=Phi_T(:,1:N_up).*exp(x2'*w*y2);
%          phi_2(:,N_up+1:N_par)=Phi_T(:,N_up+1:N_par).*exp(x2'*w*y2);
         phi_1=Phi_T;
         phi_2(:,1:N_up)=exp(x2'*w*y2)*eN2(1:N_sites,1:N_sites)*Phi_T(:,1:N_up);
         phi_2(:,N_up+1:N_par)=exp(x2'*w*y2)*eN2(N_sites+1:2*N_sites,N_sites+1:2*N_sites)*Phi_T(:,N_up+1:N_par);
         
         invO_matrix_up=inv(phi_1(:,1:N_up)'*phi_2(:,1:N_up));
         invO_matrix_dn=inv(phi_1(:,N_up+1:N_par)'*phi_2(:,N_up+1:N_par));
         
        %%  calculate the single-particle Green's function matrix for each spin:
         temp_up=phi_2(:,1:N_up)*invO_matrix_up;
         temp_dn=phi_2(:,N_up+1:N_par)*invO_matrix_dn;
         
         re_up=det(phi_1(:,1:N_up)'*phi_2(:,1:N_up));
         re_dn=det(phi_1(:,N_up+1:N_par)'*phi_2(:,N_up+1:N_par));
         re=re_up*re_dn;
         
         G_up=temp_up*phi_1(:,1:N_up)'*re;
         G_dn=temp_dn*phi_1(:,N_up+1:N_par)'*re;

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
 
end