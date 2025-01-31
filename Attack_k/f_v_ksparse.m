function [v] = f_v_ksparse(Sigma_yy,M2,Sigma_aa,sensor_candi,lambda_list) 
        n = size(Sigma_yy,2);
        C = inv(Sigma_yy + Sigma_aa);
        % D = inv(eye(n)+M2*inv(Sigma_yy))*M2;
        D = inv(eye(n)+M2*Sigma_aa)*M2;
        %% solve v
        
        a = 1/Sigma_yy(sensor_candi,sensor_candi);
        b = M2(sensor_candi,sensor_candi);
        c = C(sensor_candi,sensor_candi);
        d = D(sensor_candi,sensor_candi);
        for sensor_candi = 1:size(lambda_list,2)
            lambda = lambda_list(sensor_candi);
            inequality3 = (b - a*lambda)*(c+d) + c*d*(lambda-1); % v is positive
            inequality2 = ((a*lambda-b)*(c+d) +c*d*(1-lambda))^2 -4*(a*lambda-b)*c*d*(d-c*lambda);
            v(sensor_candi) = (inequality3 + sqrt(inequality2))./(2*c*d*(lambda*a-b));
        end
end