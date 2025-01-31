function [obj,D1,D2] = f_obj_s_onesensor(Sigma_xixi,Sigma_xiaxia,Sigma_yy,Sigma_aa,lambda,measure_ind)
    Sigma_yaya = Sigma_yy + Sigma_aa;
    n = size(Sigma_xixi,1)/2;
    m = size(Sigma_yy,1);
    Sigma_hatxhatx = Sigma_xixi(n+1:end,n+1:end);
    Sigma_hatxahatxa = Sigma_xiaxia(n+1:end,n+1:end);
    D1 = 0.5*(log(det(Sigma_hatxhatx)) - log(det(Sigma_hatxahatxa)) - n + trace(inv(Sigma_hatxhatx)*Sigma_hatxahatxa));
    D2 = 0;
    for i = 1:size(measure_ind,1)
        a = Sigma_yaya(measure_ind(i),measure_ind(i));
        b = Sigma_yy(measure_ind(i),measure_ind(i));
        D2 = D2+0.5*(log(sqrt(b)) - log(sqrt(a)) + a/(2*b) - 0.5);
    end
    obj = -D1+lambda*D2;
 
 
end