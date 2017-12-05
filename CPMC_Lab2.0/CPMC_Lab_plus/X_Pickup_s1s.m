function [E_trace,a_trace,w_trace]=X_Pickup_s1s(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,s1s_T,Phi_start)
%% invoke the main function
    %% Initialize X_RBM state, icf 2017/10/31
    X_RBM_Initialization_Pickup;
    
    N_sites=Lx*Ly*Lz;
    N_par=N_up+N_dn;

    w=rand(N_y,1);
        
    a_int=rand(N_y,N_sites);
    a=a_int;

    w_step_length=1;
    a_step_length=1;
    s1s_step_length=0.001;

    N_iteration=100;
    %
    E_trace=[];
    w_trace=[];
    w_trace(end+1,:,:) = w;
    a_trace=[];
    a_trace(end+1,:,:) = a;
    N_trace=[];
    S1S_trace=[];
    
    for j=1:N_iteration
        [E_trace,N_trace,S1S_trace,a_trace,w_trace,a,w,a_step_length,w_step_length]=X_RBM_update_Pickup_s1s(a_trace,w_trace,a,w,Phi_T,N_sites,N_y,a_step_length,w_step_length,s1s_step_length,N_up,N_dn,U,H_k,E_trace,N_trace,S1S_trace,s1s_T);
        a_step_length=a_step_length;
    end

%% Print Result
x=1:1:Lx;
y=1:1:Ly;

figure;
S1S_trace1=S1S_trace(:);
[xx,yy]=meshgrid(x,y);
zz=S1S_trace1(xx+(yy-1)*Lx);
surf(xx,yy,zz);
xlabel ('X\_site');
ylabel ('Y\_site');
zlabel ('X\_S1S');

figure;
[xx,yy]=meshgrid(x,y);
zz=N_trace(xx+(yy-1)*Lx)-N_trace(xx+(yy-1)*Lx+N_sites);
surf(xx,yy,zz);
xlabel ('X\_site');
ylabel ('Y\_site');
zlabel ('X(N\_up-N\_dn)');
% save all workplace
save ('myFile.mat');

end