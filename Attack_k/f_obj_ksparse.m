function [obj,D1,D2] = f_obj_ksparse(Sigma_xixi,Sigma_yy,Sigma_aa,lambda,A,B,K,L,C,Vd,Vn)
    Sigma_yaya = Sigma_yy + Sigma_aa;
    n = size(Sigma_xixi,1);
    m = size(Sigma_yy,1);
    % [~,~,~,~,~,Sigma_xiaxia] = cov_matrix(A,B,K,L,C,Vd,Vn,Sigma_aa);
    [~,~,~,~,~,Sigma_xiaxia] = cov_matrix(A,B,K,L,C,Vd,Vn,Sigma_aa);
    [~,~,~,~,~,Sigma_xixi] = cov_matrix(A,B,K,L,C,Vd,Vn,zeros(m,m));
 
    % D1 = 0.5*(log(det(Sigma_hatxhatx)) - log(det(Sigma_hatxahatxa)) - n + trace(inv(Sigma_hatxhatx)*Sigma_hatxahatxa));
    D1 = 0.5*(log(det(Sigma_xixi)) - log(det(Sigma_xiaxia)) - n + trace(inv(Sigma_xixi)*Sigma_xiaxia));
    D2 = 0.5*(log(det(Sigma_yy)) - log(det(Sigma_yaya)) - m + trace(inv(Sigma_yy)*Sigma_yaya));
 
    obj = - D1+lambda*D2;
 
 
end