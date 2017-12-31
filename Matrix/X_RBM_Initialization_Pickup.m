%% Trivial K State Construction
if Phi_start == 0
   H_k=H_K_pin(Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);
   [psi_nonint,E_nonint_m] = eig(H_k);
   Phi_T=horzcat(psi_nonint(:,1:N_up),psi_nonint(:,1:N_dn)); %Trivial K State.
else
   H_k=H_K_pin(Lx, Ly,Lz, kx, ky,kz, tx, ty,tz);
   Phi_T=Phi_start;
end
