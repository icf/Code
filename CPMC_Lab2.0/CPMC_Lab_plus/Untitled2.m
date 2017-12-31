% for i=1:N_L_loop
%     E_L_avg(i)=2*N_L_loop*E_L_avg(i)/(2*i);
%     E_L_err(i)=2*N_L_loop*std(E_err(i,:))/(sqrt(N_loop)*2*i);
% end
% 
% figure;
% errorbar(1./(2*(2:N_L_loop)).^2,E_L_avg(2:N_L_loop),E_L_err(2:N_L_loop)); 
% xlabel ('1/(Lx^2)');
% ylabel ('E/L');
% 
% % figure;
% % plot(1*(1:N_U_loop),E_V_ERR);
% % xlabel ('U0');
% % ylabel ('E_V_ERR');
% 
% 
% % figure;
% % plot(1*(1:N_U_loop),E_V_ERR(1:N_U_loop).*(0.5-1./sqrt(1+(8./(1:N_U_loop).^2))));
% % xlabel ('U0');
% % ylabel ('E_V_ERR*(0.5-U/sqrt(U^2+64))');
% % 
% % a=(0.5-1./sqrt(1+(8./(1:N_U_loop)).^2))
% % b=1./sqrt(1+(8./(1:N_U_loop).^2))
% 
% % save all workplace, icf 2017/9/19
% save ('myFile.mat');

% figure;
% errorbar(U,E_ave,E_err);
% xlabel ('U');
% ylabel ('E');
% 
% figure;
% errorbar(U,E_BP_ave,E_BP_err);
% xlabel ('U');
% ylabel ('E_BP');

% figure;
% plot(N_wlk,E_BP_ave-E_ED);
% xlabel ('N_wlk');
% ylabel ('E_BP_ave-E_ED');
% 
% figure;
% plot(N_wlk,E_K_BP_ave-E_K_ED_mixed);
% xlabel ('N_wlk');
% ylabel ('E_K_BP_ave-E_K_ED_mixed');
% 
% figure;
% plot(N_wlk,E_V_BP_ave-E_V_ED_mixed);
% xlabel ('N_wlk');
% ylabel ('E_V_BP_ave-E_V_ED_mixed');

% x1=bitget(uint8(3),4:-1:1)
% E_mean1=mean(E_ave1(1:7));
% E_mean2=mean(E_ave2(1:7));
% E_BP_mean1=mean(E_BP_ave1(1:7));
% E_BP_mean2=mean(E_BP_ave2(1:7));

% x=1:7;
% 
% figure;
% errorbar(x,E_ave1(x),E_err1(x),'DisplayName','E\_CPMC');
% 
% hold on
% errorbar(x,E_BP_ave1(x),E_BP_err1(x),'DisplayName','E\_CPMC\_BP');
% 
% hold on
% errorbar(x,E_ave2(x),E_err2(x),'DisplayName','E\_CPMC\_X');
% 
% hold on
% errorbar(x,E_BP_ave2(x),E_BP_err2(x),'DisplayName','E\_CPMC\_BP\_X');
% 
% ylabel ('E');
% 
% save ('myFile.mat');

Lx=16;
Ly=4;

exact_cpmc=zeros(Lx*Ly,1);

            figure;
            plot(exact_cpmc);
            xlabel ('site');
            ylabel ('hole density');








