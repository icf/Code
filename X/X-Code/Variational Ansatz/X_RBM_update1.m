function [a,w]=X_RBM_update1(a,w,Phi_T,Proj_k,N_sites,N_y,a_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k)
%
%%
flag=0;

%% Energy for any X_RBM state
E=Energy_X_RBM(a,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);

for i=1:2*N_sites+1
    delta_a=zeros(2*N_sites,1);
    if i~=2*N_sites+1
       delta_a(i,1)=a_step_length;
    end
    
    for j=1:2*N_sites+1
    for k=1:N_y+1
        delta_w=zeros(2*N_sites,N_y);
        if j~=2*N_sites+1 && k~=N_y+1
           delta_w(j,k)=w_step_length;
        end
        
        E2=Energy_X_RBM(a+delta_a,w+delta_w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
        
        if E2+E_step_length <= E  
           E=E2 
           a=a+delta_a
           w=w+delta_w
           flag=1;
           %break;
        end
    end
    end
end

if flag==0
   flag=flag
end

end