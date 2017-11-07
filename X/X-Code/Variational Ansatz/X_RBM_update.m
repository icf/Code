function [E_trace,a,w,a_step_length,w_step_length]=X_RBM_update(a,w,Phi_T,Proj_k,N_sites,N_y,a_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k,E_trace)
%
%%
flag=0;

%% Energy for any X_RBM state
[E,E_ED,E_real]=Energy_X_RBM(a,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)

for i=1:(N_up+N_dn)*0.5*N_sites+1
for i2=1:N_y+1
    delta_a=zeros((N_up+N_dn)*0.5*N_sites,N_y);
     if i~=(N_up+N_dn)*0.5*N_sites+1 && i2~=N_y+1
        delta_a(i,i2)=a_step_length;
     end
    
    for j=1:N_sites+1
    for j2=1:N_y+1
        delta_w=zeros(N_sites,N_y);
        if j~=N_sites+1 && j2~=N_y+1
           delta_w(j,j2)=w_step_length;
        end
            [E2,E_ED2,E_real2]=Energy_X_RBM(a+delta_a,w+delta_w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
        
            if E2+E_step_length <= E
               E=E2 
               E_trace(end+1) = E;
               
               E_ED=E_ED2;
               E_real=E_real2;
           
               a=a+delta_a;
               w=w+delta_w;
               flag=1;
            else
               [E2,E_ED2,E_real2]=Energy_X_RBM(a-delta_a,w-delta_w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
               if E2+E_step_length <= E
                  E=E2 
                  E_trace(end+1) = E;
                  
                  E_ED=E_ED2;
                  E_real=E_real2;
              
                  a=a-delta_a;
                  w=w-delta_w;
               end
               %break;
        end
    end
    end
end
end

if flag==0
   flag=flag
   a_step_length=0.1*a_step_length;
   w_step_length=0.1*w_step_length;
end

end