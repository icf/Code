function [E_trace,a,w,w_step_length]=X_RBM_update3_2(a,w,Phi_T,N_sites,N_y,w_step_length,E_step_length,N_up,N_dn,U,H_k,E_trace)
%
%%
flag=0;

% a_T=a';
% a_T=stblz_BP(a',N_y);
% a=a_T';
%% Energy for any X_RBM state
[E,E_ED,E_real]=Energy_X_RBM3_2(a,w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
if length(E_trace)==0
    E_trace(end+1) = E;
else
    if E<E_trace(end)
       E_trace(end+1) = E;
    end
end
    
    for j=1:N_y+1
    for j2=1:1
        delta_w=zeros(N_y,1);
        if j~=N_y+1 
           delta_w(j,j2)=w_step_length;
        end
        w0=w;
        [E22,E_ED2,E_real2]=Energy_X_RBM3_2(a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E23,E_ED2,E_real2]=Energy_X_RBM3_2(a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        %
        if E22+E_step_length <= E
           E=E22; 
           if E<E_trace(end)
              E_trace(end+1) = E;
              E=E
           end          
           w=w0+delta_w;
           flag=1;
        end
        if E23+E_step_length <= E
           E=E23; 
           if E<E_trace(end)
              E_trace(end+1) = E;
              E=E
           end          
           w=w0-delta_w;
           flag=1;
        end
        %
    end
    end


if flag==0
   if w_step_length <= 0.001 
      flag=0
      w_step_length=10000*w_step_length;
      
      a=rand(N_y,N_sites)/1000+ones(N_y,N_sites);
      w=rand(N_y,1);
      
   else
      flag=1
      w_step_length=0.1*w_step_length;
      
   end
end

end