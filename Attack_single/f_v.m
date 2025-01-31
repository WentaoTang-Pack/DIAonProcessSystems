function [v] = f_v(Sigma_yy,M2,measure_ind,lambda)
 
 
        a = 1/Sigma_yy(measure_ind,measure_ind);
        b = M2(measure_ind,measure_ind);
        v = (b^2 - lambda*a^2)/((lambda*a-b)*a*b);

     
 
end