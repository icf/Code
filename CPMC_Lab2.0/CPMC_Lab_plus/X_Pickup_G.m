function [E_trace,a_trace,w_trace]=X_Pickup_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,G_T,Phi_start)
%% invoke the main function
    %% Initialize X_RBM state, icf 2017/10/31
    X_RBM_Initialization_Pickup;
    
    N_sites=Lx*Ly*Lz;
    N_par=N_up+N_dn;

    w=zeros(N_y,1);
        
    a_int=rands(N_y,N_sites*N_sites);
    a=a_int;

    w_step_length=1;
    a_step_length=1;
    G_step_length=0.001;

    N_iteration=1000;
    %
    E_trace=[];
    w_trace=[];
    w_trace(end+1,:,:) = w;
    a_trace=[];
    a_trace(end+1,:,:) = a;
    N_trace=[];
    G_trace=[];
    
    for j=1:N_iteration
        [E_trace,N_trace,G_trace,a_trace,w_trace,a,w,a_step_length,w_step_length]=X_RBM_update_Pickup_G(a_trace,w_trace,a,w,Phi_T,N_sites,N_y,a_step_length,w_step_length,G_step_length,N_up,N_dn,U,H_k,E_trace,N_trace,G_trace,G_T);
        a_step_length=a_step_length;
    end

%% Print Result
x=1:1:Lx;
y=1:1:Ly;

figure;
G_trace1=diag(G_trace(:,1:N_sites))-diag(G_trace(:,N_sites+1:2*N_sites));
[xx,yy]=meshgrid(x,y);
zz=G_trace1(xx+(yy-1)*Lx);
surf(xx,yy,zz);
xlabel ('X\_site');
ylabel ('Y\_site');
zlabel ('X\_diag(G\_up)-diag(G\_dn)');

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