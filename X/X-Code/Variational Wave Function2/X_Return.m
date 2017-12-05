function [phi_2]=X_Return(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace)
%% Initialization
N_sites=Lx*Ly*Lz;
N_par=N_up+N_dn;
n_phi=0;

a(:,:)=a_trace(end,:,:);
w(:,1)=w_trace(end,:,1);

X_RBM_Initialization_Pickup;

%% Build up
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
         
         phi_1=Phi_T;
         phi_2(:,1:N_up,n_phi)=w(y2,1)*eN_up(:,:)*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par,n_phi)=w(y2,1)*eN_dn(:,:)*phi_1(:,N_up+1:N_par);
     end
     
end