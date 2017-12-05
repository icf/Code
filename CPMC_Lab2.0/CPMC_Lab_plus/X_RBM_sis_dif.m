function [s1s_dif] = X_RBM_sis_dif(s1s1,s1s2,N_sites)
%
%%
s1s_dif=sum(((s1s1-s1s2).*5).^2);

end