% A script to loop over multiple sets of input parameters and run a CPMC calculation for each set
%
% Huy Nguyen, Hao Shi, Jie Xu and Shiwei Zhang
% ?014 v1.0
% Package homepage: http://cpmc-lab.wm.edu
% Distributed under the <a href="matlab: web('http://cpc.cs.qub.ac.uk/licence/licence.html')">Computer Physics Communications Non-Profit Use License</a>
% Any publications resulting from either applying or building on the present package 
%   should cite the following journal article (in addition to the relevant literature on the method):
% "CPMC-Lab: A Matlab Package for Constrained Path Monte Carlo Calculations" Comput. Phys. Commun. (2014)

%% Before Loop, icf 2017/9/20
%define average loop
N_loop=10;

%define Number of L loop
N_L_loop=2;

%% system parameters:
Lx=[2:2:2*N_L_loop]; % The number of lattice sites in the x direction
Ly=1; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction

N_up=Lx./2; % The number of spin-up electrons
N_dn=Lx./2; % The number of spin-down electrons

U=4; % The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction

%% run parameters:
deltau=0.01; % The imaginary time step
N_wlk=100; % The number of random walkers
N_blksteps=200; % The number of random walk steps in each block
N_eqblk=5; % The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=10; % The number of blocks used in the measurement phase
itv_modsvd=5; % The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=10; % The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=20; % The interval between two adjacent energy measurements

%% Initialize the batch run
N_run=length(Lx); %replace argument by the parameter that needs to be looped over
E_ave=zeros(N_run,1);
E_err=zeros(N_run,1);

%% invoke the main function
for i=1:N_run
    suffix=strcat('_Lx',int2str(Lx(i))); % Set the suffix to distinguish between different runs in the same batch  
    for j=1:N_loop
        kx=2*rand()-1; % The x component of the twist angle in TABC (twist-averaging boundary condition)
        ky=0; % The y component of the twist angle in TABC
        kz=0; % The z component of the twist angle in TABC
        % call main function AFQMC_Hub
        [E_ave(i,j),E_err(i,j),savedFile]=CPMC_Lab(Lx(i),Ly,Lz,N_up(i),N_dn(i),kx,ky,kz,U,tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix); 
    end
end
%% post-run:

%% sum up, icf 2017/9/19
for i=1:N_L_loop
    E_L_avg(i)=mean(E_ave(i,:))/(2*i);
    E_L_err(i)=std(E_err(i,:))/(sqrt(N_loop)*2*i);
end

figure;
errorbar(1./(2*(1:N_L_loop)).^2,E_L_avg,E_L_err); 
xlabel ('1/(Lx^2)');
ylabel ('E/L');
% save all workplace, icf 2017/9/19
save ('myFile.mat');
%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"