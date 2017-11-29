% E_X3c3=E_trace_save(1:3);
% x=1:8;
% y=1:7;
% 
% figure;
% plot((2.*x).^2,E_Nor(x));
% % 
% hold on
% plot((2.*y).^2,E_CPMC(y));
% figure;
% plot(E_trace_save0(1:10,1));
% 
% hold on
% plot(E_trace_save1(1:10,1));
% 
% hold on
% plot(E_ave(1:10,1));
% 
% figure;
% plot((2.*y).^2,E_Nor(y,1)./E_CPMC(y,1));
% 
% hold on
% plot(E_trace_save1./E_ave);
% 

 % save all workplace
 save ('compareFile.mat');
