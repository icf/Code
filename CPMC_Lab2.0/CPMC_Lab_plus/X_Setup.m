

Lx=4;
Ly=4;
Lz=1;

kx=0;
ky=0;
kz=0;

tx=1;
ty=1;
tz=1;

N_up=7; 
N_dn=7; 

U=8;

N_y=1;  

%%
    N_sites=Lx*Ly*Lz;
    N_par=N_up+N_dn;
    H_k=H_K(Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);
    %HF
    a=0.75;
    N_it=1000;
    [Phi_start(:,:),~]=HF(Lx,Ly,Lz,N_up,N_dn,N_par,N_sites,N_it,kx,ky,kz,U,tx,ty,tz,a);  
    %
    Phi_start=Phi_start;
    %
    invO_matrix_up=inv(Phi_start(:,1:N_up)'*Phi_start(:,1:N_up));
    invO_matrix_dn=inv(Phi_start(:,N_up+1:N_par)'*Phi_start(:,N_up+1:N_par));
    temp_up=Phi_start(:,1:N_up)*invO_matrix_up;
    temp_dn=Phi_start(:,N_up+1:N_par)*invO_matrix_dn;
    G_up=temp_up*Phi_start(:,1:N_up)';
    G_dn=temp_dn*Phi_start(:,N_up+1:N_par)';
    % calculate the potential energy:
    n_int=(diag(G_up)).'*diag(G_dn);
    potentialEnergy=n_int*U;

    % calculate the kinetic energy:
    kineticEnergy=sum(sum(H_k.'.*(G_up+G_dn))); % note the element-wise multiplication

    % calculate the total energy:
    e=potentialEnergy+kineticEnergy;
    %%
    g=horzcat(G_up,G_dn);
    plot(g);
    ylabel ('HF\_G');
    %
    Phi_start=0;
    G_T=zeros(N_sites,2*N_sites);
    for i1=1:N_sites
        for i2=1:2*N_sites
        G_T(i1,i2)=g(i1,i2);
        end
    end
    [E_trace,a_trace,w_trace]=X_Pickup_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,G_T,Phi_start);
    [phi_2,w_trace]=X_Return_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace,Phi_start);
    %
    Phi_start=phi_2;
    %
    invO_matrix_up=inv(Phi_start(:,1:N_up)'*Phi_start(:,1:N_up));
    invO_matrix_dn=inv(Phi_start(:,N_up+1:N_par)'*Phi_start(:,N_up+1:N_par));
    temp_up=Phi_start(:,1:N_up)*invO_matrix_up;
    temp_dn=Phi_start(:,N_up+1:N_par)*invO_matrix_dn;
    G_up=temp_up*Phi_start(:,1:N_up)';
    G_dn=temp_dn*Phi_start(:,N_up+1:N_par)';
    %
    g=horzcat(G_up,G_dn);
    plot(g);
    ylabel ('X\_G');
    
    plot(G_T-g);
% [E_trace,a_trace,w_trace]=X_Pickup(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,phi_2);
% [phi_2,w_trace]=X_Return(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace,phi_2);
% [E_trace,a_trace,w_trace]=X_Pickup_s1s(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,s1s_T,phi_2);
% [phi_2,w_trace]=X_Return(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace,phi_2);
[E,N]=X_Check_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,w_trace,phi_2);
% pick_up=randsrc(100,1,[1:N_y; w_trace(end,:,1).^2]);



