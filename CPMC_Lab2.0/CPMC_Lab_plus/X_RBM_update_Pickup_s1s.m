function [E_trace,N_trace,s1s_trace,a_trace,w_trace,a,w,a_step_length,w_step_length]=X_RBM_update_Pickup_s1s(a_trace,w_trace,a,w,Phi_T,N_sites,N_y,a_step_length,w_step_length,s1s_step_length,N_up,N_dn,U,H_k,E_trace,N_trace,s1s_trace,s1s_T)
%
%%
flag=0;

%% Energy for any X_RBM state
[E,~,~,N,s1s]=X_RBM_Energy_X_RBM_s1s(a,w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
if length(s1s_trace)==0
    E_trace(end+1) = E;
    N_trace=N;
    s1s_trace=s1s;
    s1s_dif_st=X_RBM_sis_dif(s1s_trace,s1s_T,N_sites)
    s1s=s1s;
    s1s_dif=X_RBM_sis_dif(s1s,s1s_T,N_sites);
else
    s1s_dif=X_RBM_sis_dif(s1s,s1s_T,N_sites);
    s1s_dif_st=X_RBM_sis_dif(s1s_trace,s1s_T,N_sites);
    if s1s_dif<s1s_dif_st
       E_trace(end+1) = E;
       a_trace(end+1,:,:) = a;
       w_trace(end+1,:,:) = w;
       N_trace=N;
       s1s_trace=s1s;
       s1s_dif_st=X_RBM_sis_diff(s1s_trace,s1s_T,N_sites);
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
        [E21,~,~,N21,s1s21]=X_RBM_Energy_X_RBM_s1s(a+delta_a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E22,~,~,N22,s1s22]=X_RBM_Energy_X_RBM_s1s(a-delta_a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E23,~,~,N23,s1s23]=X_RBM_Energy_X_RBM_s1s(a+delta_a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E24,~,~,N24,s1s24]=X_RBM_Energy_X_RBM_s1s(a-delta_a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        %
        s1s_dif21=X_RBM_sis_dif(s1s21,s1s_T,N_sites);
        if s1s_dif21+s1s_step_length <= s1s_dif
           s1s=s1s21; 
           E=E21;
           s1s_dif=X_RBM_sis_dif(s1s,s1s_T,N_sites);
           if s1s_dif<s1s_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0+delta_a;
              w_trace(end+1,:,:) = w0+delta_w;
              N_trace=N21;
              s1s_trace=s1s21;
              s1s_dif_st=X_RBM_sis_dif(s1s_trace,s1s_T,N_sites)
              E=E
           end          
           a=a0+delta_a;
           w=w0+delta_w;
           flag=1;
        end
        s1s_dif22=X_RBM_sis_dif(s1s22,s1s_T,N_sites);
        if s1s_dif22+s1s_step_length <= s1s_dif
           s1s=s1s22; 
           E=E22;
           s1s_dif=X_RBM_sis_dif(s1s,s1s_T,N_sites);
           if s1s_dif<s1s_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0+delta_a;
              w_trace(end+1,:,:) = w0+delta_w;
              N_trace=N22;
              s1s_trace=s1s22;
              s1s_dif_st=X_RBM_sis_dif(s1s_trace,s1s_T,N_sites)
              E=E
           end          
           a=a0-delta_a;
           w=w0+delta_w;
           flag=1;
        end
        s1s_dif23=X_RBM_sis_dif(s1s23,s1s_T,N_sites);
        if s1s_dif23+s1s_step_length <= s1s_dif
           s1s=s1s23; 
           E=E23;
           s1s_dif=X_RBM_sis_dif(s1s,s1s_T,N_sites);
           if s1s_dif<s1s_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0+delta_a;
              w_trace(end+1,:,:) = w0+delta_w;
              N_trace=N23;
              s1s_trace=s1s23;
              s1s_dif_st=X_RBM_sis_dif(s1s_trace,s1s_T,N_sites)
              E=E
           end          
           a=a0+delta_a;
           w=w0-delta_w;
           flag=1;
        end
        s1s_dif23=X_RBM_sis_dif(s1s23,s1s_T,N_sites);
        if s1s_dif23+s1s_step_length <= s1s_dif
           s1s=s1s24; 
           E=E24;
           s1s_dif=X_RBM_sis_dif(s1s,s1s_T,N_sites);
           if s1s_dif<s1s_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0+delta_a;
              w_trace(end+1,:,:) = w0+delta_w;
              N_trace=N24;
              s1s_trace=s1s24;
              s1s_dif_st=X_RBM_sis_dif(s1s_trace,s1s_T,N_sites)
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
      flag=0;
      a_step_length=10000*a_step_length;
      w_step_length=10000*w_step_length;
      
      a=rand(N_y,N_sites);
      w=rand(N_y,1);
   else
      flag=1;
      a_step_length=0.1*a_step_length;
      w_step_length=0.1*w_step_length;
   end
end

end