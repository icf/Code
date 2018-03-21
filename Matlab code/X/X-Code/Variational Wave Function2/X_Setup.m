

Lx=4;
Ly=4;
Lz=1;

kx=0;
ky=0;
kz=0;

tx=1;
ty=1;
tz=1;

N_up=Ly; 
N_dn=Ly; 

U=4;

N_y=4;  

%%
[E_trace,a_trace,w_trace]=X_Pickup(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y);
[phi_2]=X_Return(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz,N_up,N_dn,U,N_y,E_trace,a_trace,w_trace);

