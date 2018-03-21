%% started X_RBM state


%% Trivial K State Construction
H_k=H_K(Lx(i), Ly(i),Lz, kx, ky,kz, tx, ty,tz);
[psi_nonint,E_nonint_m] = eig(H_k);
Phi_T=horzcat(psi_nonint(:,1:N_up(i)),psi_nonint(:,1:N_dn(i))); %Trivial K State.

Proj_k = expm(-0.5*u*H_k);
