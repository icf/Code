function [E_trace,N_trace,a,w,a_step_length,w_step_length]=X_RBM_update3_2(a,w,Phi_T,N_sites,N_y,a_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k,E_trace,N_trace)
%
%%
flag=0;

% a_T=a';
% a_T=stblz_X(a',N_y);
% a=a_T';
%% Energy for any X_RBM state
[E,E_ED,E_real,N]=X_RBM_Energy_X_RBM3_2(a,w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
if length(E_trace)==0
    E_trace(end+1) = E;
    N_trace=N;
else
    if E<E_trace(end)
       E_trace(end+1) = E;
       N_trace=N;
    end
end

for i=1:N_y+1
for i2=1:N_sites+1    
    delta_a=zeros(N_y,N_sites);
     if i~=N_y+1 && i2~=N_sites+1 
        delta_a(i,i2)=a_step_length;
     end
    
    for j=1:N_y+1
    for j2=1:1
        delta_w=zeros(N_y,1);
        if j~=N_y+1 
           delta_w(j,j2)=w_step_length;
        end
        a0=a;
        w0=w;
        [E21,E_ED2,E_real2,N21]=X_RBM_Energy_X_RBM3_2(a+delta_a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E22,E_ED2,E_real2,N22]=X_RBM_Energy_X_RBM3_2(a-delta_a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E23,E_ED2,E_real2,N23]=X_RBM_Energy_X_RBM3_2(a+delta_a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E24,E_ED2,E_real2,N24]=X_RBM_Energy_X_RBM3_2(a-delta_a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        %
        if E21+E_step_length <= E
           E=E21; 
           if E<E_trace(end)
              E_trace(end+1) = E;
              N_trace=N21;
              E=E
           end          
           a=a0+delta_a;
           w=w0+delta_w;
           flag=1;
        end
        if E22+E_step_length <= E
           E=E22; 
           if E<E_trace(end)
              E_trace(end+1) = E;
              N_trace=N22;
              E=E
           end          
           a=a0-delta_a;
           w=w0+delta_w;
           flag=1;
        end
        if E23+E_step_length <= E
           E=E23; 
           if E<E_trace(end)
              E_trace(end+1) = E;
              N_trace=N23;
              E=E
           end          
           a=a0+delta_a;
           w=w0-delta_w;
           flag=1;
        end
        if E24+E_step_length <= E
           E=E24; 
           if E<E_trace(end)
              E_trace(end+1) = E;
              N_trace=N24;
              E=E
           end          
           a=a0-delta_a;
           w=w0-delta_w;
           flag=1;
        end
        %
    end
    end
end
end


if flag==0
   if a_step_length <= 0.001 
      flag=0
      a_step_length=10000*a_step_length;
      w_step_length=10000*w_step_length;
      
      a=rand(N_y,N_sites);
      w=rand(N_y,1);
   else
      flag=1
      a_step_length=0.1*a_step_length;
      w_step_length=0.1*w_step_length;
   end
end

end