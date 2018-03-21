

Lx=16;
Ly=4;
Lz=1;

kx=0;
ky=0;
kz=0;

tx=1;
ty=1;
tz=1;

N_up=28; 
N_dn=28; 

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
    E1=potentialEnergy+kineticEnergy;
    %%
    g=horzcat(G_up,G_dn);
    figure;
    plot(g);
    ylabel ('HF\_g');
    %
    Phi_start=0;
    G_T=zeros(N_sites,2*N_sites);
    for i1=1:N_sites
        for i2=1:2*N_sites
        G_T(i1,i2)=g(i1,i2);
        end
    end
    %% Matrix_X, icf, 2017/12/16 
    X_RBM_Initialization_Pickup;
    %
    [Phi_chang_up, Phi_chang_dn]=Matrix_X(G_T, N_up, N_dn, N_par, N_sites);
    phi_2(:,1:N_up)=Phi_chang_up;
    phi_2(:,N_up+1:N_par)=Phi_chang_dn;
    %
    Phi_start=phi_2;
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
    E2=potentialEnergy+kineticEnergy;
    %
    g=horzcat(G_up,G_dn);
    g=g;
    figure;
    plot(g);
    ylabel ('Matrix\_X\_G');

    Ggerr=sum(sum(((G_T-g).*100).^2))



