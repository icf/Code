% 
% E_Nor2x16=E_trace_save(1:7);
%     
x=1:10;
y=1:7;
% z=1:5;
% k=1:5;
% j=1:3;
% i=1:2;
% % 
figure;
plot(2.*x,E_Nor2x(x)./CPMC_X2x(x),'DisplayName','E\_Nor2x/E\_CPMC');

hold on
plot(2.*x,E_Nor2x8(x)./CPMC_X2x(x),'DisplayName','E\_Nor2x8/E\_CPMC');

hold on
plot(2.*x,E_X2x4(x)./CPMC_X2x(x),'DisplayName','E\_X4/E\_CPMC'); 

hold on
plot(2.*x,E_X2x8(x)./CPMC_X2x(x),'DisplayName','E\_X8/E\_CPMC');

hold on
plot(2.*y,E_X2x16(y)./CPMC_X2x(y),'DisplayName','E\_X16/E\_CPMC');

% figure;
% plot(2.*x,(E_X2x8(x)-E_X2x4(x))./CPMC_X2x(x),'DisplayName','E_X2x8-E_X2x4');
% 
% hold on
% plot((2.*i).^2,E_X16(i)./E_CPMC(i),'DisplayName','E\_X16/E\_CPMC');
%     
% hold on
% plot((2.*y).^2,E_CPMC(y)./E_CPMC(y),'DisplayName','E\_CPMC/E\_CPMC'); 
% 
% hold on
% plot((2.*z).^2,E_X2c2(z)./E_CPMC(z),'DisplayName','E\_X2c2/E\_CPMC'); 
% 
% hold on
% plot((2.*k).^2,E_X3c3(k)./E_CPMC(k),'DisplayName','E\_X3c3/E\_CPMC'); 

%     
%     figure;
%     plot(abs(E_ave(3:25,1)-E_trace_save4(3:25,1))./abs(E_trace_save3(3:25,1)-E_trace_save4(3:25,1)));
%     hold on
%     plot(:);
%     
%     figure;
%     plot(E_X4./E_CPMC);
%     
%     hold on
%     plot(E_Nor(1:7,1)./E_CPMC);

% save all workplace
save ('compareFile.mat');