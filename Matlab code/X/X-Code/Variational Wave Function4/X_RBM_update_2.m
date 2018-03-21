function [E_trace,a,b,w,a_step_length,b_step_length,w_step_length]=X_RBM_update_2(a,b,w,Phi_T,Proj_k,N_sites,N_y,a_step_length,b_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k,E_trace)
%
%%
flag=0;

%% Energy for any X_RBM state
[E,E_ED,E_real]=Energy_X_RBM_2(a,b,w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
if length(E_trace)==0
    E_trace(end+1) = E;
else
    if E<E_trace(end)
       E_trace(end+1) = E;
    end
end


for i=1:N_y+1
for i2=1:(N_up+N_dn)*0.5*N_sites+1
    delta_a=zeros(N_y,(N_up+N_dn)*0.5*N_sites);
     if i~=N_y+1 && i2~=(N_up+N_dn)*0.5*N_sites+1 
        delta_a(i,i2)=a_step_length;
     end
    
    for j=1:N_y+1
    for j2=1:1
        delta_w=zeros(N_y,1);
        if j~=N_y+1 
           delta_w(j,j2)=w_step_length;
        end
        for k=1:1:(N_up+N_dn)*0.5*N_sites+1
            delta_b=zeros((N_up+N_dn)*0.5*N_sites,1);
            if k~=(N_up+N_dn)*0.5*N_sites+1 
               delta_b(k,1)=b_step_length;
            end
            [E2,E_ED2,E_real2]=Energy_X_RBM_2(a+delta_a,b+delta_b,w+delta_w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
        
            if E2+E_step_length <= E
               E=E2; 
               if E<E_trace(end)
                  E_trace(end+1) = E;
                  E=E
               end
               
               E_ED=E_ED2;
               E_real=E_real2;
           
               a=a+delta_a;
               w=w+delta_w;
               b=b+delta_b;
               flag=1;
            else
               [E2,E_ED2,E_real2]=Energy_X_RBM_2(a-delta_a,b+delta_b,w-delta_w,Phi_T,Proj_k,N_sites,N_y,N_up,N_dn,U,H_k);
               if E2+E_step_length <= E
                  E=E2; 
                  if E<E_trace(end)
                     E_trace(end+1) = E;
                     E=E
                  end
                  
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
   if a_step_length <= 0.001 && w_step_length <= 0.001
      flag=1;

      a_step_length=10000*a_step_length;
      w_step_length=10000*w_step_length;
      b_step_length=10000*b_step_length;
      
      a=rand(N_y,(N_up+N_dn)*0.5*N_sites);
      b=ones((N_up+N_dn)*0.5*N_sites,1);
      w=rand(N_y,1);
   else
      flag=1;
      a_step_length=0.1*a_step_length;
      w_step_length=0.1*w_step_length;
      b_step_length=0.1*b_step_length;
   end

end



end