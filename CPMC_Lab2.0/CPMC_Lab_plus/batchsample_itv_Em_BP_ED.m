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
Lx=2; % The number of lattice sites in the x direction
Ly=2; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction

N_up=1; % The number of spin-up electrons
N_dn=1; % The number of spin-down electrons

kx=0; % The x component of the twist angle in TABC (twist-averaging boundary condition)
ky=0; % The y component of the twist angle in TABC
kz=0; % The z component of the twist angle in TABC

U=4; % The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction

%% run parameters:
deltau=0.01; % The imaginary time step
N_wlk=[100:50:100]; % The number of random walkers
N_blksteps=300; % The number of random walk steps in each block
N_eqblk=5; % The number of blocks used to equilibrate the random walk before energy measurement takes place
N_blk=10; % The number of blocks used in the measurement phase
itv_modsvd=5; % The interval between two adjacent modified Gram-Schmidt re-orthonormalization of the random walkers. No re-orthonormalization if itv_modsvd > N_blksteps
itv_pc=5; % The interval between two adjacent population controls. No population control if itv_pc > N_blksteps
itv_Em=20; % The interval between two adjacent energy measurements

%% Initialize the batch run
N_run=length(N_wlk); %replace argument by the parameter that needs to be looped over
E_ave=zeros(N_run,1);
E_err=zeros(N_run,1);
E_BP_ave=zeros(N_run,1);
E_BP_err=zeros(N_run,1);
E_K_BP_ave=zeros(N_run,1);
E_K_BP_err=zeros(N_run,1);
E_V_BP_ave=zeros(N_run,1);
E_V_BP_err=zeros(N_run,1);
E_ED=zeros(N_run,1);
E_K_ED_mixed=zeros(N_run,1);
E_V_ED_mixed=zeros(N_run,1);

%% invoke the main function
for i=1:N_run
    suffix=strcat('_itv_Em',int2str(N_wlk(i))); % Set the suffix to distinguish between different runs in the same batch    
    % call main function AFQMC_Hub
    [E_ave(i),E_err(i),E_BP_ave(i),E_BP_err(i),E_V_BP_ave(i),E_V_BP_err(i),E_K_BP_ave(i),E_K_BP_err(i),E_ED(i),E_V_ED_mixed(i),E_K_ED_mixed(i),savedFile]=CPMC_Lab_ED(Lx,Ly,Lz,N_up,N_dn,kx,ky,kz,U,tx,ty,tz,deltau,N_wlk(i),N_blksteps,N_eqblk,N_blk,itv_modsvd,itv_pc,itv_Em,suffix);    
end
%% post-run:
% plot energy vs different run parameters
figure;
errorbar(N_wlk,E_ave,E_err);
xlabel ('N_wlk');
ylabel ('E');

hold on
errorbar(N_wlk,E_BP_ave,E_BP_err);
xlabel ('N_wlk');
ylabel ('E_BP');

hold on
plot(N_wlk,E_ED);
xlabel ('N_wlk');
ylabel ('E_ED');

figure;
errorbar(N_wlk,abs(E_BP_ave-E_ED),E_BP_err);
xlabel ('N\_wlk');
ylabel ('E\_BP\_err');
%

figure;
errorbar(N_wlk,E_K_BP_ave,E_K_BP_err);
xlabel ('N_wlk');
ylabel ('E_K_BP');

hold on
plot(N_wlk,E_K_ED_mixed);
xlabel ('N_wlk');
ylabel ('E_K_ED_mixed');
%

figure;
errorbar(N_wlk,E_V_BP_ave,E_V_BP_err);
xlabel ('N_wlk');
ylabel ('E_V_BP');

hold on
plot(N_wlk,E_V_ED_mixed);
xlabel ('N_wlk');
ylabel ('E_V_ED_mixed');
%

figure;
errorbar(N_wlk,abs(E_K_BP_ave-E_K_ED_mixed),E_K_BP_err);
xlabel ('N_wlk');
ylabel ('E\_K\_BP\_err');
%

figure;
errorbar(N_wlk,abs(E_V_BP_ave-E_V_ED_mixed),E_V_BP_err);
xlabel ('N_wlk');
ylabel ('E\_V\_BP\_err');


% save all workplace, icf 2017/9/19
save ('myFile.mat');
%% Explanation of saved quantities:
% E: the array of energy of each block
% time: The total computational time
% E_nonint_v: the non-interacting energy levels of the system
% Phi_T: the trial wave function
% For other saved quantities, type "help CPMC_Lab"