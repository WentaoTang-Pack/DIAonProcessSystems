function [bestsensor] = f_bestSensor_ksparse(sensor_candi,Sigma_yy,M2,Sigma_aa,Sigma_xixi,A,B,K,L,C,Vd,Vn,lambda_condition)
    n = size(Sigma_yy,2);
    for i = 1:size(sensor_candi,2)
        [lambda_min,lambda_max] = f_lambda_lim_ksparse(Sigma_yy,M2,Sigma_aa,sensor_candi(i));
        if lambda_min > lambda_condition || lambda_max < lambda_condition 
            continue
        end
        disp(i);
        v(i) = f_v_ksparse(Sigma_yy,M2,Sigma_aa,sensor_candi(i),lambda_condition);   
        temp = zeros(n,n); temp(sensor_candi(i),sensor_candi(i)) = v(i); % v is for different lambda
        Sigma_aa_temp = Sigma_aa + temp;
        [obj_i_lambda(i),~,~] = f_obj_ksparse(Sigma_xixi,Sigma_yy,Sigma_aa_temp,lambda_condition,A,B,K,L,C,Vd,Vn);
    end
    obj_i_lambda_temp = obj_i_lambda(obj_i_lambda ~= 0);
        [val,~] = min(obj_i_lambda_temp);  
        bestsensor = find(obj_i_lambda == val);
 
end