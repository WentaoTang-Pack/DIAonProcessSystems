function [bestsensor] = f_bestSensor_ksparse_20250113(sensor_candi,Sigma_yy,M2,Sigma_aa,Sigma_xixi,A,B,K,L,C,Vd,Vn,lambda)
    n = size(Sigma_yy,2);
    for i = 1:size(sensor_candi,2)
        % lambda_list = lambda_min+0.01:0.1:lambda_max - 0.01;
        v = f_v_ksparse(Sigma_yy,M2,Sigma_aa,sensor_candi(i),lambda);
 
            temp = zeros(n,n); temp(sensor_candi(i),sensor_candi(i)) = v; % v is for different lambda
            Sigma_aa_temp = Sigma_aa + temp;
            [obj_i_sensor(i),~,~] = f_obj_ksparse(Sigma_xixi,Sigma_yy,Sigma_aa_temp,lambda,A,B,K,L,C,Vd,Vn);
    end
    [~,index] = min(obj_i_sensor);  
 
    bestsensor = sensor_candi(index);
 
 
end