function [Phi,n_old_out]=HF(Lx,Ly,Lz,N_up,N_dn,N_par,N_sites,N_it,kx,ky,kz,U,tx,ty,tz,a)
%Mean Field (UHF approximation), icf 2017/10/16
%% Initialization
H_k=zeros(N_sites,N_sites);
psi_nonint=zeros(N_sites,N_sites);
E_nonint_m=zeros(N_sites,N_sites);
Phi_T=zeros(N_sites,N_par);
n_new_in=zeros(2*N_sites);
n_old_in=zeros(2*N_sites);
n_old_out=zeros(2*N_sites);
err=zeros(N_it);
%counter
N_it_counter=1:1:N_it;
%plot 
x=1:1:Lx;
y=1:1:Ly;

H_k=H_K_pin(Lx,Ly,Lz,kx,ky,kz,tx,ty,tz);
[psi_nonint,~] = eig(H_k);
Phi_T=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn));

n_new_in=HF_n(Phi_T,N_sites,N_up,N_dn,N_par);
% n_new_in=[1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,1,0,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0];
%n_new_in=[0.25,0.25,0.25,0.25,0.25,0.25,0.25,0.25];
%% HF iteration
for i=1:N_it
    n_old_in=n_new_in;
%     Phi=HF_H(n_old_in,H_k,N_sites,N_up,N_dn,U);
    Phi=HF_H_pin(n_old_in,H_k,N_sites,N_up,N_dn,U,Lx,Ly,Lz,kx,ky,kz,tx,ty,tz);
    n_old_out=HF_n(Phi,N_sites,N_up,N_dn,N_par);
    n_new_in=(1-a)*n_old_in+a*n_old_out;
    err(i)=sum((n_old_out-n_old_in).^2);
    if N_it-i<1
       figure;
       [xx,yy]=meshgrid(x,y);
       zz=n_old_out(xx+(yy-1)*Lx)-n_old_out(xx+(yy-1)*Lx+N_sites);
       surf(xx,yy,zz);
       xlabel ('X\_site');
       ylabel ('Y\_site');
       zlabel ('n\_up-n\_dn');
       
       figure;
       [xx,yy]=meshgrid(x,y);
       zz=n_old_out(xx+(yy-1)*Lx)+n_old_out(xx+(yy-1)*Lx+N_sites);
       surf(xx,yy,zz);
       xlabel ('X\_site');
       ylabel ('Y\_site');
       zlabel ('n\_up+n\_dn');
    end
end

figure;
plot(N_it_counter,err);
xlabel ('n\_it');
ylabel ('err');

end