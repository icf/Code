function [a,w]=X_RBM_update4(a,w,b,Phi_T,Proj_k,N_sites,N_y,a_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k)
%
%%
flag=0;

%% Energy for any X_RBM state
[E,E_ED,E_real]=Energy_X_RBM4(a,w,b,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k)

for i=1:N_sites+1
for i2=1:1
    delta_a=zeros(N_sites,1);
     if i~=N_sites+1 
        delta_a(i,i2)=a_step_length;
     end
    
    for j=1:N_sites+1
    for j2=1:N_y+1
        delta_w=zeros(N_sites,N_y);
        if j~=N_sites+1 && j2~=N_y+1
           delta_w(j,j2)=w_step_length;
        end
        for k=1:2*N_sites+1
            delta_b=zeros(2*N_sites,1);
            if k~=2*N_sites+1 
               delta_b(k,1)=w_step_length;
            end
            [E2,E_ED2,E_real2]=Energy_X_RBM4(a+delta_a,w+delta_w,b+delta_b,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
        
            if E2+E_step_length <= E
               E=E2 
               E_ED=E_ED2;
               E_real=E_real2;
           
               a=a+delta_a;
               w=w+delta_w;
               b=b+delta_b;
               flag=1;
            else
               [E2,E_ED2,E_real2]=Energy_X_RBM4(a-delta_a,w-delta_w,b+delta_b,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
               if E2+E_step_length <= E
                  E=E2 
                  E_ED=E_ED2;
                  E_real=E_real2;
              
                  a=a-delta_a;
                  w=w-delta_w;
                  b=b-delta_b;
               end
               %break;
            end
        end
    end
    end
end
end

if flag==0
   flag=flag
end

end