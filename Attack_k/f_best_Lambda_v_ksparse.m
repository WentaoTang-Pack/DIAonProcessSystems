function [opt_lambda,opt_v] = f_best_Lambda_v_ksparse(Sigma_yy,M2,Sigma_aa,ind,Sigma_xixi,A,B,K,L,C,Vd,Vn)
    n = size(Sigma_yy,2);
    [lambda_min,lambda_max] = f_lambda_lim_ksparse(Sigma_yy,M2,Sigma_aa,ind);
    lambda_list = lambda_min+0.01:0.1:lambda_max - 0.01;
    v = f_v_ksparse(Sigma_yy,M2,Sigma_aa,ind,lambda_list);
    for j = 1:size(lambda_list,2)
        lambda = lambda_list(j);
        temp = zeros(n,n); temp(ind,ind) = v(j); % v is for different lambda
        Sigma_aa_temp = Sigma_aa + temp;
        [obj_i(j),~,~] = f_obj_ksparse(Sigma_xixi,Sigma_yy,Sigma_aa_temp,lambda,A,B,K,L,C,Vd,Vn);
    end
    [~,index] = min(obj_i); obj_i = [];
    opt_lambda = lambda_list(index);
    opt_v = v(index);
 
 
end