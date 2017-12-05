function [phi_2,w_trace]=X_Return_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace,Phi_start)
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
                 eN_up(i,i)=1;
                 eN_dn(i,i)=1;
         end  
         for i=1:N_sites
             for j=1:N_sites
                 eN_up(i,j)=eN_up(i,j)+w(y2,1)*exp(temp_T(i+(j-1)*N_sites,1));
                 eN_dn(i,j)=eN_dn(i,j)+w(y2,1)*exp(-1*temp_T(i+(j-1)*N_sites,1));
             end
         end  
         eN_up=X_RBM_stblz_X(eN_up,N_sites);
         eN_up=X_RBM_stblz_X(eN_dn,N_sites);

         phi_1=Phi_T;
         phi_2(:,1:N_up,n_phi)=eN_up(:,:)*phi_1(:,1:N_up);
         phi_2(:,N_up+1:N_par,n_phi)=eN_dn(:,:)*phi_1(:,N_up+1:N_par);
         
         w(y2,1)=w(y2,1)/w_sum;
     end
     w_trace(end,:,1)=w(:,1);
     
end