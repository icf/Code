function [Phi_chang_up1, Phi_chang_dn1, prob1_up, prob1_dn, Phi_chang_up2, Phi_chang_dn2, prob2_up, prob2_dn]=Matrix_X_mixed(G_T, N_up, N_dn, N_par, N_sites, Lx, Ly, G2_T_up, G2_T_dn)
%% initialization
    I=zeros(N_sites,N_sites);
    for i=1:N_sites
        I(i,i)=1;
    end
%
    a=0;
%   
    G_T_up=G_T(:,1:N_sites);
    G_T_dn=G_T(:,1+N_sites:2*N_sites);

    G_T_up=((1-a).*G_T_up+a.*G2_T_up);
    G_T_up=(G_T_up+G_T_up')./2;
    G_T_dn=((1-a).*G_T_dn+a.*G2_T_dn);
    G_T_dn=(G_T_dn+G_T_dn')./2;
    
%% PROB1
    [psi_nonint_up,E_nonint_m] = eig(G_T_up);
    E_nonint_v_up=diag(E_nonint_m)
    E_pick=[];
    order=1;
    prob1_up=0;
    for i=1:N_sites
        if abs(E_nonint_v_up(i)) >= 0.5 && order <= N_up
            order=order+1;
            E_pick(end+1)=i;
            prob1_up=prob1_up+abs(E_nonint_v_up(i));
        end
    end
    prob1_up=prob1_up/N_up
    E_pick=E_pick
    Phi_2_up=horzcat(psi_nonint_up(:,E_pick));
    %
    [psi_nonint_dn,E_nonint_m] = eig(G_T_dn);
    E_nonint_v_dn=diag(E_nonint_m)
    E_pick=[];
    order=1;
    prob1_dn=0;
    for i=1:N_sites
        if abs(E_nonint_v_dn(i)) >= 0.5 && order <= N_dn
            order=order+1;
            E_pick(end+1)=i;
            prob1_dn=prob1_dn+abs(E_nonint_v_dn(i));
        end
    end
    prob1_dn=prob1_dn/N_dn
    E_pick=E_pick
    Phi_2_dn=horzcat(psi_nonint_dn(:,E_pick));

Phi_chang_up1=Phi_2_up;
Phi_chang_dn1=Phi_2_dn;

%% PROB2
    [psi_nonint_up,E_nonint_m] = eig(G_T_up);
    E_nonint_v_up=diag(E_nonint_m)
    E_pick=[];
    order=1;
    prob2_up=0;
    for i=1:N_sites
        if abs(E_nonint_v_up(i)) <= 0.5 && i>2 &&  order <= N_up
            order=order+1;
            E_pick(end+1)=i;
            prob2_up=prob2_up+abs(E_nonint_v_up(i));
        end
    end
    prob2_up=prob2_up/N_up
    E_pick=E_pick
    Phi_2_up=horzcat(psi_nonint_up(:,E_pick));
    %
    [psi_nonint_dn,E_nonint_m] = eig(G_T_dn);
    E_nonint_v_dn=diag(E_nonint_m)
    E_pick=[];
    order=1;
    prob2_dn=0;
    for i=1:N_sites
        if abs(E_nonint_v_dn(i)) <= 0.5 && i>2 && order <= N_dn
            order=order+1;
            E_pick(end+1)=i;
            prob2_dn=prob2_dn+abs(E_nonint_v_dn(i));
        end
    end
    prob2_dn=prob2_dn/N_dn
    E_pick=E_pick
    Phi_2_dn=horzcat(psi_nonint_dn(:,E_pick));
%% POST RUN
Phi_chang_up2=Phi_2_up;
Phi_chang_dn2=Phi_2_dn;

prob1_up=prob1_up/(prob1_up+prob2_up);
prob2_up=prob2_up/(prob1_up+prob2_up);

prob1_dn=prob1_dn/(prob1_dn+prob2_dn);
prob2_dn=prob2_dn/(prob1_dn+prob2_dn);

G_up=Phi_chang_up1*Phi_chang_up1'*prob1_up+Phi_chang_up2*Phi_chang_up2'*prob2_up;
Err_Phi_2_up=sum(sum(  ((G_up-G_T_up).*100).^2 ))

G_dn=Phi_chang_dn1*Phi_chang_dn1'*prob1_dn+Phi_chang_dn2*Phi_chang_dn2'*prob2_dn;
Err_Phi_2_dn=sum(sum(  ((G_dn-G_T_dn).*100).^2 ))

end