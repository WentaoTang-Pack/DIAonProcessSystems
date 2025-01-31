function [lambda_min,lambda_max] = f_lambda_lim_single(Sigma_yy,M2,measure_ind)
 
            a = 1/Sigma_yy(measure_ind,measure_ind);
            b = M2(measure_ind,measure_ind);
            lambda_min = b/a;
            lambda_max = b^2/a^2;
 
         
 
 
end