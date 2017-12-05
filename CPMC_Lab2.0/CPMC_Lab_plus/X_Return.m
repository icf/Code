function [phi_2,w_trace]=X_Return(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace,Phi_start)
%% Initialization
N_sites=Lx*Ly*Lz;
N_par=N_up+N_dn;
n_phi=0;
w_sum=0;

a(:,:)=a_trace(end,:,:);
w(:,1)=w_trace(end,:,1);

X_RBM_Initialization_Pickup;

%% Build up
for y1=1:N_y
    w_sum=w_sum+(w(y1,1))^2;
end
w_sum=sqrt(w_sum);

     for y2=1:N_y %binary order
         n_phi=n_phi+1;
         eN_up=zeros(N_sites,N_sites);
         eN_dn=zeros(N_sites,N_sites);
         
         if N_y==1
            temp(:,1)=a(:);
            temp_T=temp;
         else
            temp=a(y2,:);
            temp_T=temp';
         end
         
         for i=1:N_sites
             eN_up(i,i)=exp(temp_T(i,1));
             eN_dn(i,i)=exp(-1*temp_T(i,1));
         end  
         
         phi_1=Phi_T;
         phi_2(:,1:N_up,n_phi)=eN_up(:,:)*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par,n_phi)=eN_dn(:,:)*phi_1(:,N_up+1:N_par);
         
         phi_2(:,1:N_up,n_phi)=X_RBM_stblz_X(phi_2(:,1:N_up,n_phi),N_up);
         phi_2(:,N_up+1:N_par,n_phi)=X_RBM_stblz_X(phi_2(:,N_up+1:N_par,n_phi),N_dn);
         
         w(y2,1)=w(y2,1)/w_sum;
     end
     w_trace(end,:,1)=w(:,1);
     
end