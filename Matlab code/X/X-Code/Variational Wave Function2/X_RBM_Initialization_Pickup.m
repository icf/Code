%% Trivial K State Construction
H_k=H_K(Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);
[psi_nonint,E_nonint_m] = eig(H_k);
Phi_T=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn)); %Trivial K State.
