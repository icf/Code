function [Phi_chang_up, Phi_chang_dn]=Matrix_X(G_T, N_up, N_dn, N_par, N_sites)
%% initialization
    I=zeros(N_sites,N_sites);
    for i=1:N_sites
        I(i,i)=1;
    end

    G_T_up=G_T(:,1:N_sites);
    G_T_dn=G_T(:,1+N_sites:2*N_sites);
    G1_up=rands(N_sites,N_sites);
    G1_up=G1_up+G1_up';
    G1_dn=rands(N_sites,N_sites);
    G1_dn=G1_dn+G1_dn';
%     for i=1:N_sites
%         G1_up(i,i)=1;
%         G1_dn(i,i)=1;
%     end
    %
    it_run=1000;
%%
    [psi_nonint,E_nonint_m] = eig(G_T_up);
    E_nonint_v=diag(E_nonint_m)
    E_pick=[];
    order=1;
    for i=1:N_sites
        if abs(E_nonint_v(i)) >= 0.9 && order <= N_up
            order=order+1;
            E_pick(end+1)=i;
        end
    end
    E_pick=E_pick
    % assemble the non-interacting single-particle orbitals into a Slater determinant:
    Phi_2_up=horzcat(psi_nonint(:,E_pick));
    invO_matrix_up=inv(Phi_2_up'*Phi_2_up);
    Gg_up=Phi_2_up*Phi_2_up';
    Err_Phi_2_up=sum(sum(  ((Gg_up-G_T_up).*100).^2 ))
    %
    [psi_nonint,E_nonint_m] = eig(G_T_dn);
    E_nonint_v=diag(E_nonint_m)
    E_pick=[];
    order=1;
    for i=1:N_sites
        if abs(E_nonint_v(i)) >= 0.9 && order <= N_dn
            order=order+1;
            E_pick(end+1)=i;
        end
    end
    % assemble the non-interacting single-particle orbitals into a Slater determinant:
    E_pick=E_pick
    Phi_2_dn=horzcat(psi_nonint(:,E_pick));
    invO_matrix_up=inv(Phi_2_dn'*Phi_2_dn);
    Gg_dn=Phi_2_dn*Phi_2_dn';
    Err_Phi_2_up=sum(sum(  ((Gg_dn-G_T_dn).*100).^2 ))

%%
% for i=1:it_run
% %     invO_matrix_up1=inv(Phi_T_up'*G1_up'*G1_up*Phi_T_up);
% %     invO_matrix_dn1=inv(Phi_T_dn'*G1_dn'*G1_dn*Phi_T_dn);
%     invO_matrix_up=inv(Phi_T_up*Phi_T_up'+I);
%     invO_matrix_dn=inv(Phi_T_dn*Phi_T_dn'+I);
%     inv_G1_up_T=G1_up;
%     inv_G1_dn_T=G1_dn;
%     %
%     G2_up=(G_T_up+I)*G1_up*invO_matrix_up;
%     G2_dn=(G_T_dn+I)*G1_dn*invO_matrix_dn;
%     
%     G2_up=X_RBM_stblz_X(G2_up,N_sites);
%     G2_dn=X_RBM_stblz_X(G2_dn,N_sites);
%     
%     G1_up=G2_up+0.75*G1_up;
%     G1_dn=G2_dn+0.75*G1_dn;
%     
%     Err_G_up(i)=sum(sum(abs(  G2_up-(G_T_up+I)*G2_up*invO_matrix_up )));
% 
% end
% 
% figure;
% plot(Err_G_up);

Phi_chang_up=Phi_2_up;
Phi_chang_dn=Phi_2_dn;

end