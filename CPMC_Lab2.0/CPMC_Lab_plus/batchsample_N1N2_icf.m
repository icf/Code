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
%define deltaU and Number of deltaU
N_deltaU=1; % N_deltaU has to be 1 or E_K_ERR won't give the right estimation.
deltaU=0.08;
N_loop=10;

%define Number of U loop
N_U_loop=8;

%% run parameters:
deltau=0.01; % The imaginary time step
N_wlk=100; % The number of random walkers
N_blksteps=200; % The number of random walk steps in each block
N_eqblk=5; % The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=10; % The number of blocks used in the measurement phase
itv_modsvd=5; % The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=10; % The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=20; % The interval between two adjacent energy measurements

%% system parameters:
Lx=2; % The number of lattice sites in the x direction
Ly=1; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction

N_up=1; % The number of spin-up electrons
N_dn=1; % The number of spin-down electrons

kx=0; % The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=0; % The y component of the twist angle in TABC
kz=0; % The z component of the twist angle in TABC

%% Start K Loop, icf 2017/9/20
for i_U=1:N_U_loop
U0=1*i_U;
display(U0);
%% Start Loop, icf 2017/9/19
for i_N_loop=1:N_loop
%% kx 
U=U0; % The on-site repulsion strength in the Hubbard Hamiltonian
U12=[U0-N_deltaU*deltaU:deltaU:U0+N_deltaU*deltaU];
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction

%% Initialize the batch run
N_run=length(U12); %replace argument by the parameter that needs to be looped over
E_ave=zeros(N_run,1);
E_err=zeros(N_run,1);

%% invoke the main function
for i=1:N_run
    suffix=strcat('_U12',int2str(U12(i))); % Set the suffix to distinguish between different runs in the same batch    
    % call main function AFQMC_Hub
    [E_ave(i),E_err(i),savedFile]=CPMC_Lab_N1N2(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,U12(i),tx,ty,tz,deltau,N_wlk,N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix);    
end
%% post-run:
% plot energy vs different run parameters
% figure;
% errorbar(tx,E_ave,E_err);
% xlabel ('tx');
% ylabel ('E');
% calculate and plot E_NiNi, icf 2017/9/20
for i=1:N_deltaU
    E_N1N2(i,i_N_loop)=(E_ave(2+2*N_deltaU-i)-E_ave(i))./(2*(1+N_deltaU-i)*deltaU);
end

%% end loop
end
% plot E_NiNi 
for i=1:N_deltaU
    E_N1N2_avg(i_U,i)=mean(E_N1N2(i,:));
    E_N1N2_err(i_U,i)=std(E_N1N2(i,:))/sqrt(N_loop);
end
E_N1N2_AVG(i_U)=mean(E_N1N2_avg(i_U,:));
E_N1N2_ERR(i_U)=mean(E_N1N2_err(i_U,:));
%% end U loop, icf 2017/9/20
end

figure;
errorbar(1*(1:N_U_loop),E_N1N2_AVG,E_N1N2_ERR); 
xlabel ('U0');
ylabel ('E_N1N2');
% save all workplace, icf 2017/9/19
save ('myFile.mat');
%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"