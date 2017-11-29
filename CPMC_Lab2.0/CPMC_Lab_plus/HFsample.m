%% system parameters:
Lx=16; % The number of lattice sites in the x direction
Ly=16; % The number of lattice sites in the y direction
Lz=1; % The number of lattice sites in the z direction !only plot 2D Density Picture, so Lz=1.

N_up=; % The number of spin-up electrons
N_dn=540; % The number of spin-down electrons

% kx=2*rand(1)-1; % The x component of the twist angle in TABC (twist-averaging boundary condition)
% ky=kx*Ly/Lx;
kx=0;
ky=0;
%ky=2*rand(1)-1; % The y component of the twist angle in TABC
kz=0; % The z component of the twist angle in TABC

U=9.0; % The on-site repulsion strength in the Hubbard Hamiltonian
tx=1; % The hopping amplitude between nearest-neighbor sites in the x direction
ty=1; % The hopping amplitude between nearest neighbor sites in the y direction
tz=1; % The hopping amplitude between nearest neighbor sites in the z direction

N_it=800;
a=0.75;

%% Initialize the batch run
N_sites=Lx*Ly*Lz;
N_par=N_up+N_dn;
N_run=length(kx); %replace argument by the parameter that needs to be looped over
Phi=zeros(N_sites,N_par,N_run);
n=zeros(2*N_sites,N_run);
N=zeros(2*N_sites);
n_ave=zeros(2*N_sites);
%plot 
x=1:1:Lx;
y=1:1:Ly;

%% invoke the main function
for i=1:N_run
    % call main function HF
    [Phi(:,:,i),n(:,i)]=HF(Lx,Ly,Lz,N_up,N_dn,N_par,N_sites,N_it,kx(i),ky(i),kz,U,tx,ty,tz,a);    
end
%% post-run:
% save all workplace, icf 2017/9/19
save ('HFmyFile.mat');
% plot energy vs different run parameters
for i=1:N_run
    N=N+n(:,i);
end
n_ave=N./N_run;

figure;
[xx,yy]=meshgrid(x,y);
zz=n_ave(xx+(yy-1)*Lx)+n_ave(xx+(yy-1)*Lx+N_sites);
surf(xx,yy,zz);
xlabel ('X\_site');
ylabel ('Y\_site');
zlabel ('n+\_ave');

% figure;
% [xx,yy]=meshgrid(x,y);
% zz=n_ave(xx+(yy-1)*Lx);
% surf(xx,yy,zz);
% xlabel ('X\_site');
% ylabel ('Y\_site');
% zlabel ('n\_up\_ave');
% 
% figure;
% [xx,yy]=meshgrid(x,y);
% zz=n_ave(xx+(yy-1)*Lx+N_sites);
% surf(xx,yy,zz);
% xlabel ('X\_site');
% ylabel ('Y\_site');
% zlabel ('n\_dn\_ave');