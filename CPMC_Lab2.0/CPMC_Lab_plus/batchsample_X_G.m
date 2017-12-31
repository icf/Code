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
deltau=0.01; % The imaginary time step
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
E_K_BP_ave=zeros(N_run,1);
E_K_BP_err=zeros(N_run,1);
E_V_BP_ave=zeros(N_run,1);
E_V_BP_err=zeros(N_run,1);

%% invoke the main function
for i=1:N_run
    suffix=strcat('_itv_Em',int2str(N_wlk(i))); % Set the suffix to distinguish between different runs in the same batch    
    % call main function AFQMC_Hub
    [E_ave1(i),E_err1(i),E_BP_ave1(i),E_BP_err1(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),G_BP_ave1(i,:,:),G_BP_err1(i,:,:),savedFile]=CPMC_Lab_G(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix);    
%
    N_sites=Lx*Ly*Lz;
    figure;
    G_BP_ave_2D(:,1:N_sites)=G_BP_ave1(i,:,1:N_sites);
    G_BP_ave_2D(:,N_sites+1:2*N_sites)=G_BP_ave1(i,:,N_sites+1:2*N_sites);
    G_BP_ave_p=diag(G_BP_ave_2D(:,1:N_sites))-diag(G_BP_ave_2D(:,N_sites+1:2*N_sites));
    [xx,yy]=meshgrid(x,y);
    zz=G_BP_ave_p(xx+(yy-1)*Lx);
    surf(xx,yy,zz);
    xlabel ('X\_site');
    ylabel ('Y\_site');
    zlabel ('BP\_diag(G\_up)-diag(G\_dn)');

    %% X_s1s
    N_y=1
    %
    N_sites=Lx*Ly*Lz;
    %
    Phi_start=0;
    G_T=zeros(N_sites,2*N_sites);
    for i1=1:N_sites
        for i2=1:2*N_sites
        G_T(i1,i2)=G_BP_ave1(i,i1,i2);
        end
    end
    %
    [E_trace,a_trace,w_trace]=X_Pickup_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,G_T,Phi_start);
    [phi_2,w_trace]=X_Return_G(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace,Phi_start);  
    %
    [E_ave2(i),E_err2(i),E_BP_ave2(i),E_BP_err2(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),savedFile]=CPMC_Lab_BP_X(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix,phi_2);  
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