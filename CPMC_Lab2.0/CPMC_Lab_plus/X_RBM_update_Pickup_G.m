function [E_trace,N_trace,G_trace,a_trace,w_trace,a,w,a_step_length,w_step_length]=X_RBM_update_Pickup_G(a_trace,w_trace,a,w,Phi_T,N_sites,N_y,a_step_length,w_step_length,G_step_length,N_up,N_dn,U,H_k,E_trace,N_trace,G_trace,G_T)
%
%%
flag=0;

%% Energy for any X_RBM state
[E,~,~,N,G]=X_RBM_Energy_X_RBM_G(a,w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
if length(G_trace)==0
    E_trace(end+1) = E
    N_trace=N;
    G_trace=G;
    G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites)
    G=G;
    G_dif=X_RBM_G_dif(G,G_T,N_sites);
else
    G_dif=X_RBM_G_dif(G,G_T,N_sites);
    G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites);
    if G_dif<G_dif_st
       E_trace(end+1) = E;
       a_trace(end+1,:,:) = a;
       w_trace(end+1,:,:) = w;
       N_trace=N;
       G_trace=G;
       G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites)
    end
end

for i=1:N_y
for i2=1:N_sites*N_sites+1    
    delta_a=zeros(N_y,N_sites*N_sites);
     if i2~=N_sites*N_sites+1 
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
        [E21,~,~,N21,G21]=X_RBM_Energy_X_RBM_G(a+delta_a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E22,~,~,N22,G22]=X_RBM_Energy_X_RBM_G(a-delta_a,w+delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E23,~,~,N23,G23]=X_RBM_Energy_X_RBM_G(a+delta_a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        [E24,~,~,N24,G24]=X_RBM_Energy_X_RBM_G(a-delta_a,w-delta_w,Phi_T,N_sites,N_y,N_up,N_dn,U,H_k);
        %
        G_dif21=X_RBM_G_dif(G21,G_T,N_sites);
        if G_dif21+G_step_length <= G_dif
           G=G21; 
           E=E21;
           G_dif=X_RBM_G_dif(G,G_T,N_sites);
           if G_dif<G_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0+delta_a;
              w_trace(end+1,:,:) = w0+delta_w;
              N_trace=N21;
              G_trace=G21;
              G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites)
              E=E
           end          
           a=a0+delta_a;
           w=w0+delta_w;
           flag=1;
        end
        G_dif22=X_RBM_G_dif(G22,G_T,N_sites);
        if G_dif22+G_step_length <= G_dif
           G=G22; 
           E=E22;
           G_dif=X_RBM_G_dif(G,G_T,N_sites);
           if G_dif<G_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0-delta_a;
              w_trace(end+1,:,:) = w0+delta_w;
              N_trace=N22;
              G_trace=G22;
              G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites)
              E=E
           end          
           a=a0-delta_a;
           w=w0+delta_w;
           flag=1;
        end
        G_dif23=X_RBM_G_dif(G23,G_T,N_sites);
        if G_dif23+G_step_length <= G_dif
           G=G23; 
           E=E23;
           G_dif=X_RBM_G_dif(G,G_T,N_sites);
           if G_dif<G_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0+delta_a;
              w_trace(end+1,:,:) = w0-delta_w;
              N_trace=N23;
              G_trace=G23;
              G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites)
              E=E
           end          
           a=a0+delta_a;
           w=w0-delta_w;
           flag=1;
        end
        G_dif23=X_RBM_G_dif(G23,G_T,N_sites);
        if G_dif23+G_step_length <= G_dif
           G=G24; 
           E=E24;
           G_dif=X_RBM_G_dif(G,G_T,N_sites);
           if G_dif<G_dif_st
              E_trace(end+1) = E;
              a_trace(end+1,:,:) = a0-delta_a;
              w_trace(end+1,:,:) = w0-delta_w;
              N_trace=N24;
              G_trace=G24;
              G_dif_st=X_RBM_G_dif(G_trace,G_T,N_sites)
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
   if a_step_length <= 0.01 
      flag=0;
      a_step_length=1000*a_step_length;
      w_step_length=1000*w_step_length;
      
%       a=rands(N_y,N_sites*N_sites);
      w=rands(N_y,1);
   else
      flag=1;
      a_step_length=0.1*a_step_length;
      w_step_length=0.1*w_step_length;
   end
end

end