function [M2] = f_M2(A,B,K,L,C,Sigma_xx,Sigma_hatxhatx,Sigma_xhatx,Sigma_hatxx)
    M2 = zeros(size(Sigma_xx));
    max = 1000;
    X = inv(Sigma_hatxhatx - Sigma_hatxx*inv(Sigma_xx)*Sigma_xhatx);
    for k = 0:max
        M2 = M2 + ((A-L*C + B*K)')^k*X*(A-L*C+B*K)^k;
    end
    M2 = L'*M2*L;
    % S = M2 - lambda*M1;
    % P = Sigma_yy*M2;
    % Q = -lambda*M1*inv(M2);
end