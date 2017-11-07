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
a=ones((N_up+N_dn)*0.5*N_sites,N_y);
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
