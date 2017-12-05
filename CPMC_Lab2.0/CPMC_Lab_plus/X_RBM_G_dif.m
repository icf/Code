function [G_dif] = X_RBM_G_dif(G1,G2,N_sites)
%
%%
G_dif=sum(sum( ((G1-G2).*100).^2 ));

end