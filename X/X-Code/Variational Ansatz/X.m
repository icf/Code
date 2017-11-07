%% system parameters:
Lx=6;
Ly=1;
Lz=1;

kx=0;
ky=0;
kz=0;

tx=1;
ty=1;
tz=1;

N_up=1; 
N_dn=1; 
N_par=N_up+N_dn;

u=-1;                 % Time Evolution.
U=4;

%% running parameters:
N_y=1;                % Number of y (Hidden Field).
N_sites=Lx*Ly*Lz;
%a=[0.781000000000000;0.771000000000000;0.780000000000000;0.792000000000000;0.812000000000000;0.781000000000000;0.991000000000000;1.02000000000000;1;1.02000000000000;1.05000000000000;1;0.990000000000000;1.04000000000000;1.05000000000000;1.03000000000000;1.02000000000000;1.03000000000000];
a=ones((N_up+N_dn)*0.5*N_sites,N_y);
%w=[1.86100000000000;-0.0110000000000000;0.00300000000000000;-0.0150000000000000;0.00300000000000000;-0.0100000000000000];
w=ones(N_sites,N_y);

w_step_length=0.1;
a_step_length=0.1;
E_step_length=0.0001;

N_iteration=1000;

%% Initialize X_RBM state, icf 2017/10/31
X_RBM_Initialization;

%% Minimize Energy, icf 2017/10/31
% Get Energy and update for any X_RBM state
E_trace=[];
for i=1:N_iteration
    [E_trace,a,w,a_step_length,w_step_length]=X_RBM_update(a,w,Phi_T,Proj_k,N_sites,N_y,a_step_length,w_step_length,E_step_length,N_up,N_dn,U,H_k,E_trace);
    figure;
    plot(abs(E_trace));
end

%% Print Result
% save all workplace
save ('myFile.mat');
