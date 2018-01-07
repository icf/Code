% A script to loop over multiple sets of input parameters and run a CPMC calculation for each set
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ?014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% system parameters:
Lx=4; % The number of lattice sites in the x direction
Ly=4; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction
%
N_sites=Lx*Ly*Lz;

N_up=7; % The number of spin-up electrons
N_dn=7; % The number of spin-down electrons

kx=0; % The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=0; % The y component of the twist angle in TABC
kz=0; % The z component of the twist angle in TABC

U=8; % The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction

%% run parameters:
sc_run=0;
%
deltau=0.05; % The imaginary time step
N_wlk=[100:1:100]; % The number of random walkers
N_blksteps=100; % The number of random walk steps in each block
N_eqblk=5; % The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=10; % The number of blocks used in the measurement phase
itv_modsvd=5; % The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=10; % The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=20; % The interval between two adjacent energy measurements

%
x=1:1:Lx;
y=1:1:Ly;

%% Initialize the batch run
N_run=length(N_wlk); %replace argument by the parameter that needs to be looped over
E_ave=zeros(N_run,1);
E_err=zeros(N_run,1);
E_BP_ave=zeros(N_run,1);
E_BP_err=zeros(N_run,1);
S1S_BP_ave=zeros(N_run,Lx*Ly*Lz);
S1S_BP_err=zeros(N_run,Lx*Ly*Lz);
%
E_ave1=zeros(N_run,1);
E_err1=zeros(N_run,1);
E_BP_ave1=zeros(N_run,1);
E_BP_err1=zeros(N_run,1);
S1S_BP_ave1=zeros(N_run,Lx*Ly*Lz);
S1S_BP_err1=zeros(N_run,Lx*Ly*Lz);
%
E_ave2=zeros(N_run,1);
E_err2=zeros(N_run,1);
E_BP_ave2=zeros(N_run,1);
E_BP_err2=zeros(N_run,1);
S1S_BP_ave2=zeros(N_run,Lx*Ly*Lz);
S1S_BP_err2=zeros(N_run,Lx*Ly*Lz);
%
E_ave3=zeros(N_run,1);
E_err3=zeros(N_run,1);
E_BP_ave3=zeros(N_run,1);
E_BP_err3=zeros(N_run,1);
S1S_BP_ave3=zeros(N_run,Lx*Ly*Lz);
S1S_BP_err3=zeros(N_run,Lx*Ly*Lz);
%
E_ave4=zeros(N_run,1);
E_BP_ave4=zeros(N_run,1);
%
E_K_BP_ave=zeros(N_run,1);
E_K_BP_err=zeros(N_run,1);
E_V_BP_ave=zeros(N_run,1);
E_V_BP_err=zeros(N_run,1);

%
H_k=H_K_pin(Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);

%% invoke the main function
for i=1:N_run
    suffix=strcat('_itv_Em',int2str(N_wlk(i))); % Set the suffix to distinguish between different runs in the same batch    
    % call main function AFQMC_Hub
    [E_ave(i),E_err(i),E_BP_ave(i),E_BP_err(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),G_BP_ave(i,:,:),G_BP_err(i,:,:),savedFile]=CPMC_Lab_G_pin(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix);    

   %% post run
    G_BP_ave_up_temp(:,:)=G_BP_ave(i,:,1:N_sites);
    Density_BP_ave_1D_up(1:N_sites)=diag(G_BP_ave_up_temp);
    G_BP_ave_dn_temp(:,:)=G_BP_ave(i,:,N_sites+1:2*N_sites);
    Density_BP_ave_1D_dn(1:N_sites)=diag(G_BP_ave_dn_temp);

    Density_Ly_ave=0;
    for j=1:Ly
        Density_Ly_ave=Density_Ly_ave+1.-(Density_BP_ave_1D_up(1+Lx*(j-1):Lx+Lx*(j-1))+Density_BP_ave_1D_dn(1+Lx*(j-1):Lx+Lx*(j-1)));
        figure;
        plot(1.-(Density_BP_ave_1D_up(1+Lx*(j-1):Lx+Lx*(j-1))+Density_BP_ave_1D_dn(1+Lx*(j-1):Lx+Lx*(j-1))));
        xlabel ('site');
        ylabel ('hole density');
        
    end
    Density_Ly_ave=Density_Ly_ave/Ly;
    figure;
    plot(Density_Ly_ave);
    xlabel ('site');
    ylabel ('hole density ave');
    
    %% S.C. Green, icf 2017/12/7
    N_y=1
    %
    N_sites=Lx*Ly*Lz;
    N_par=N_up+N_dn;
    %HF
    a=0.75;
    N_it=1000;
    [Phi_start(:,:),~]=HF(Lx,Ly,Lz,N_up,N_dn,N_par,N_sites,N_it,kx,ky,kz,U,tx,ty,tz,a);  
    %
    G_T=zeros(N_sites,2*N_sites);
    for i1=1:N_sites
        for i2=1:2*N_sites
        G_T(i1,i2)=G_BP_ave(i,i1,i2);
        end
    end
    %    
    G2_T_up=G_T(:,1:N_sites);
    G2_T_dn=G_T(:,N_sites+1:2*N_sites);

    %
    for j=1:sc_run
        [Phi_chang_up1, Phi_chang_dn1, prob1_up, prob1_dn, Phi_chang_up2, Phi_chang_dn2, prob2_up, prob2_dn]=Matrix_X_mixed( G_T, N_up, N_dn, N_par, N_sites, Lx, Ly, G2_T_up, G2_T_dn);
        %1
        phi_2(:,1:N_up)=Phi_chang_up1;
        phi_2(:,N_up+1:N_par)=Phi_chang_dn1;
        [E_ave1(i),E_err1(i),E_BP_ave1(i),E_BP_err1(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),G_BP_ave1(i,:,:),G_BP_err1(i,:,:),savedFile]=CPMC_Lab_G_X_pin(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix,phi_2);  
        G_T1=zeros(N_sites,2*N_sites);
        for i1=1:N_sites
            for i2=1:2*N_sites
                G_T1(i1,i2)=G_BP_ave1(i,i1,i2);
            end
        end
        %2
        phi_2(:,1:N_up)=Phi_chang_up2;
        phi_2(:,N_up+1:N_par)=Phi_chang_dn2;
        [E_ave2(i),E_err2(i),E_BP_ave2(i),E_BP_err2(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),G_BP_ave2(i,:,:),G_BP_err2(i,:,:),savedFile]=CPMC_Lab_G_X_pin(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix,phi_2);  
        G_T2=zeros(N_sites,2*N_sites);
        for i1=1:N_sites
            for i2=1:2*N_sites
                G_T2(i1,i2)=G_BP_ave2(i,i1,i2);
            end
        end
        %3
        G_T(:,1:N_sites)=prob1_up*G_T1(:,1:N_sites)+prob2_up*G_T2(:,1:N_sites);
        G_T(:,1+N_sites:2*N_sites)=prob1_dn*G_T1(:,1+N_sites:2*N_sites)+prob2_dn*G_T2(:,1+N_sites:2*N_sites);
        [Phi_chang_up3, Phi_chang_dn3, ~, ~, ~, ~, ~, ~]=Matrix_X_mixed( G_T, N_up, N_dn, N_par, N_sites, Lx, Ly, G2_T_up, G2_T_dn);
        phi_2(:,1:N_up)=Phi_chang_up3;
        phi_2(:,N_up+1:N_par)=Phi_chang_dn3;
        [E_ave3(i),E_err3(i),E_BP_ave3(i),E_BP_err3(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),G_BP_ave3(i,:,:),G_BP_err3(i,:,:),savedFile]=CPMC_Lab_G_X_pin(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix,phi_2);  
        G_T=zeros(N_sites,2*N_sites);
        for i1=1:N_sites
            for i2=1:2*N_sites
                G_T(i1,i2)=G_BP_ave3(i,i1,i2);
            end
        end
        %
        G_BP_ave_up_temp(:,:)=G_BP_ave3(i,:,1:N_sites);
        Density_BP_ave_1D_up(1:N_sites)=diag(G_BP_ave_up_temp);
        G_BP_ave_dn_temp(:,:)=G_BP_ave3(i,:,N_sites+1:2*N_sites);
        Density_BP_ave_1D_dn(1:N_sites)=diag(G_BP_ave_dn_temp);

        Density_Ly_ave=0;
        for j1=1:Ly
            Density_Ly_ave=Density_Ly_ave+1.-(Density_BP_ave_1D_up(1+Lx*(j1-1):Lx+Lx*(j1-1))+Density_BP_ave_1D_dn(1+Lx*(j1-1):Lx+Lx*(j1-1)));
            figure;
            plot(1.-(Density_BP_ave_1D_up(1+Lx*(j1-1):Lx+Lx*(j1-1))+Density_BP_ave_1D_dn(1+Lx*(j1-1):Lx+Lx*(j1-1))));
            xlabel ('site');
            ylabel (['hole density_X, sc\_run',num2str(j)]);   
        end
    
        Density_Ly_ave=Density_Ly_ave/Ly;
        figure;
        plot(Density_Ly_ave);
        xlabel ('site');
        ylabel (['hole density ave_X, sc\_run',num2str(j)]);
        %
        E_ave4=prob1_up*E_ave1+prob2_up*E_ave2;
        E_BP_ave4=prob1_up*E_BP_ave1+prob2_up*E_BP_ave2;
        E_sc_save(j)=E_ave4;
        E_sc_BP_save(j)=E_BP_ave4;
    end
    
end

%% post-run:
    

%% save all workplace, icf 2017/9/19
save ('myFile.mat');
%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"