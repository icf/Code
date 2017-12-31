function [Phi_chang_up, Phi_chang_dn]=Matrix_X(G_T, N_up, N_dn, N_par, N_sites, Lx, Ly, G2_T_up, G2_T_dn)
%% initialization
    I=zeros(N_sites,N_sites);
    for i=1:N_sites
        I(i,i)=1;
    end
%
    a=0.5;
%   

% G_T(:,1:N_sites)=(G_T(:,1:N_sites)+G_T(:,1:N_sites)')./2;
% G_T(:,1+N_sites:2*N_sites)=(G_T(:,1+N_sites:2*N_sites)+G_T(:,1+N_sites:2*N_sites)')./2;
    temp_up_diag=diag(G_T(:,1:N_sites));
    temp_dn_diag=diag(G_T(:,1+N_sites:2*N_sites));
    G_T_up_diag=zeros(N_sites,N_sites);
    G_T_dn_diag=zeros(N_sites,N_sites);
    for i=1:N_sites
        G_T_up_diag(i,i)=temp_up_diag(i);
        G_T_dn_diag(i,i)=temp_dn_diag(i);
    end

    G_T_up=((1-a).*G_T_up_diag+a.*G2_T_up);
    G_T_up=(G_T_up+G_T_up')./2;
    G_T_dn=((1-a).*G_T_dn_diag+a.*G2_T_dn);
    G_T_dn=(G_T_dn+G_T_dn')./2;
    
    it_run=1000;
    
%%
    [psi_nonint_up,E_nonint_m] = eig(G_T_up);
    E_nonint_v_up=diag(E_nonint_m)
    E_pick=[];
    order=1;
    for i=1:N_sites
        if abs(E_nonint_v_up(i)) >= 0.5 && order <= N_up
            order=order+1;
            E_pick(end+1)=i;
        end
    end
    E_pick=E_pick
    Phi_2_up=horzcat(psi_nonint_up(:,E_pick));
%     invO_matrix_up=inv(Phi_2_up'*Phi_2_up);
    %
    [psi_nonint_dn,E_nonint_m] = eig(G_T_dn);
    E_nonint_v_dn=diag(E_nonint_m)
    E_pick=[];
    order=1;
    for i=1:N_sites
        if abs(E_nonint_v_dn(i)) >= 0.5 && order <= N_dn
            order=order+1;
            E_pick(end+1)=i;
        end
    end
    E_pick=E_pick
    Phi_2_dn=horzcat(psi_nonint_dn(:,E_pick));
%     invO_matrix_up=inv(Phi_2_dn'*Phi_2_dn);

%% Matrix_variation
%     Phi_2_chang_up=zeros(N_sites,N_up);
%     for i=1:N_up
%         Phi_2_chang_up(i,i)=1;
%     end
%     Phi_2_chang_dn=zeros(N_sites,N_dn);
%     for i=1:N_dn
%         Phi_2_chang_dn(i,i)=1;
%     end
%     chang_up=zeros(N_sites,N_sites);
%     chang_dn=zeros(N_sites,N_sites);
%     
%     for i=1:N_sites
%         chang_up(i,i)=E_nonint_v_up(i);
%         chang_dn(i,i)=E_nonint_v_dn(i);
%     end
%     
% %
% Err=3;
% for i=1:it_run
%     for j=1:N_sites
%         for k=1:N_up
%         Phi_2_chang_up2=Phi_2_chang_up;
%         Phi_2_chang_up2(j,k)=Phi_2_chang_up2(j,k)+1;
%         Phi_2_chang_up2=X_RBM_stblz_X(Phi_2_chang_up2,N_up);
%         
%         Err2=sum(sum(abs((Phi_2_chang_up2*Phi_2_chang_up2'-chang_up))));
%         if Err2 < Err
%            Phi_2_chang_up=Phi_2_chang_up2;
%            Err=Err2
%         end
%         end
%     end
% end
% 
% Err=3;
% for i=1:it_run
%     for j=1:N_sites
%         for k=1:N_dn
%         Phi_2_chang_dn2=Phi_2_chang_dn;
%         Phi_2_chang_dn2(j,k)=Phi_2_chang_dn2(j,k)+1;
%         Phi_2_chang_dn2=X_RBM_stblz_X(Phi_2_chang_dn2,N_up);
%         
%          Err2=sum(sum(abs((Phi_2_chang_dn2*Phi_2_chang_dn2'-chang_dn))));
%         if Err2 < Err
%            Phi_2_chang_dn=Phi_2_chang_dn2;
%            Err=Err2
%         end
%         end
%     end
% end
% %
%     
%     temp_up=diag(Phi_2_chang_up*Phi_2_chang_up');
%     chang_up_Err=sum(abs(temp_up-E_nonint_v_up))
%     temp_dn=diag(Phi_2_chang_dn*Phi_2_chang_dn');
%     chang_dn_Err=sum(abs(temp_dn-E_nonint_v_dn))
%     
%     Phi_2_up=psi_nonint_up*Phi_2_chang_up;
%     Phi_2_dn=psi_nonint_dn*Phi_2_chang_dn;
    
%% plot
%         Density_up(1:N_sites)=diag(G_up);
%         Density_dn(1:N_sites)=diag(G_dn);
% 
%         Density_Ly_ave=0;
%         
%         for j1=1:Ly
%             Density_Ly_ave=Density_Ly_ave+1.-(Density_up(1+Lx*(j1-1):Lx+Lx*(j1-1))+Density_dn(1+Lx*(j1-1):Lx+Lx*(j1-1)));
%             figure;
%             plot(1.-(Density_up(1+Lx*(j1-1):Lx+Lx*(j1-1))+Density_dn(1+Lx*(j1-1):Lx+Lx*(j1-1))));
%             xlabel ('site');
%             ylabel (['hole density_X, Matrix']);   
%         end
%     
%         Density_Ly_ave=Density_Ly_ave/Ly;
%         figure;
%         plot(Density_Ly_ave);
%         xlabel ('site');
%         ylabel (['hole density ave_X, Matrix']);

    G_up=Phi_2_up*Phi_2_up';
    Err_Phi_2_up=sum(sum(  ((G_up-G_T_up).*100).^2 ))

    G_dn=Phi_2_dn*Phi_2_dn';
    Err_Phi_2_dn=sum(sum(  ((G_dn-G_T_dn).*100).^2 ))

Phi_chang_up=Phi_2_up;
Phi_chang_dn=Phi_2_dn;

end