function [E_trace,a_trace,w_trace]=X_Pickup(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,Phi_start)
%
%% invoke the main function
    %% Initialize X_RBM state, icf 2017/10/31
    X_RBM_Initialization_Pickup;
    
    N_sites=Lx*Ly*Lz;
    N_par=N_up+N_dn;

    w=rand(N_y,1);
        
    a_int=rand(N_y,N_sites);
    a=a_int;

    w_step_length=0.1;
    a_step_length=1;
    E_step_length=0.001;

    N_iteration=25;
    %
    E_trace=[];
    w_trace=[];
    w_trace(end+1,:,:) = w;
    a_trace=[];
    a_trace(end+1,:,:) = a;
    N_trace=[];
    
    for j=1:N_iteration
        [E_trace,N_trace,a_trace,w_trace,a,w,a_step_length,w_step_length]=X_RBM_update_Pickup(a_trace,w_trace,a,w,Phi_T,N_sites,N_y,a_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k,E_trace,N_trace);
        a_step_length=a_step_length;
    end

%% Print Result
% save all workplace
save ('myFile.mat');

end