function [Sigma_xx, Sigma_xhatx,Sigma_hatxx,Sigma_hatxhatx,Sigma_yy,Sigma_xixi] = cov_matrix(A,B,K,L,C,Vd,Vn,Sigma_aa)
    n = size(A,1);
    %% do it in onego
    Block = [A,B*K;L*C, A-L*C + B*K];
    Q = [Vd, zeros(n,n); zeros(n,n), L*(Vn + Sigma_aa)*L'];
    Sigma_xixi = dlyap(Block, Q);
    Sigma_xx = Sigma_xixi(1:n,1:n);
    Sigma_hatxx = Sigma_xixi(n+1:2*n,1:n);
    Sigma_xhatx = Sigma_hatxx';
    Sigma_hatxhatx = Sigma_xixi(n+1:2*n,n+1:2*n);
    %% Sigmam_yy
    Sigma_yy = C*Sigma_xx*C' + Vn;
    
end