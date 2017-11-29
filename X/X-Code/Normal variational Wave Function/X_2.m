%% system parameters:
Lx=2;
Ly=[2:2:20];
Lz=1;

kx=0;
ky=0;
kz=0;

tx=1;
ty=1;
tz=1;

N_up=Lx.*Ly./2; 
N_dn=Lx.*Ly./2;

U=4;

%% running parameters:
N_y=4;                % Number of y (Hidden Field).

%% Minimize Energy, icf 2017/10/31
% Get Energy and update for any X_RBM state
N_run=length(Ly); %replace argument by the parameter that needs to be looped over
E_trace_save=zeros(N_run,1);

%% invoke the main function
for i=1:N_run
   
    %% Initialize X_RBM state, icf 2017/10/31
    X_RBM_Initialization_2;
    
    N_sites=Lx*Ly(i)*Lz;
    N_par=N_up(i)+N_dn(i);

    a=rand(N_y,N_sites)/1000+ones(N_y,N_sites);

    w=rand(N_y,1);

    w_step_length=1;
    E_step_length=0.001;

    N_iteration=25;
    %
    E_trace=[];
    N_trace=[];
    for j=1:N_iteration
        [E_trace,a,w,w_step_length]=X_RBM_update3_2(a,w,Phi_T,N_sites,N_y,w_step_length,E_step_length,N_up(i),N_dn(i),U,H_k,E_trace);
        w_step_length=w_step_length
    end
    
    figure;
    plot(E_trace);
    
    E_trace_save(i,1)=E_trace(end);

end

x=1:length(Ly);

figure;
plot((2.*x).^2,E_trace_save(x,1));


%% Print Result
% save all workplace
save ('myFile.mat');