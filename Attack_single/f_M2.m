function [M2] = f_M2(A,B,K,L,C,Sigma_xx,Sigma_hatxhatx,Sigma_xhatx,Sigma_hatxx,Sigma_xixi)
    M2 = zeros(size(Sigma_xixi));
    max = 1000;
    X = inv(Sigma_xixi);
    BigMatrix = [A B*K; L*C A-L*C+B*K];
    for k = 0:max
        M2 = M2 + (BigMatrix')^k*X*(BigMatrix)^k;
    end
    M2 = [zeros(size(Sigma_xx)) L']*M2*[zeros(size(Sigma_xx));L];
    % S = M2 - lambda*M1;
    % P = Sigma_yy*M2;
    % Q = -lambda*M1*inv(M2);
end