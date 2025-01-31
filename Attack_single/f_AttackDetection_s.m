function [P_d,P_f] = f_AttackDetection_s(Sigma_yy,v,measure_ind)
    tau = 1;
    N = 200000;
    a = Sigma_yy(measure_ind,measure_ind)+v;
    b = Sigma_yy(measure_ind,measure_ind);
    ya = normrnd(0, a, N,1);
    y = normrnd(0, b, N,1);
 
    % Test on probability of detection
    H1 = 0;
    H0 = 0;
    for i = 1:N  
        L1 = exp(-0.5*ya(i)^2/(2*a))./ sqrt(a*2*pi);
        L2 = exp(-0.5*ya(i)^2/(2*b))./ sqrt(b*2*pi);
        if L1/L2 > tau
            H1 = H1+1;
        else
            H0 = H0+1;
        end
    end
    P_d = H1/N;
    % Test on false almarm
    H1 = 0;
    H0 = 0;
    for i = 1:N %% use compromised data to test
        L1 = exp(-0.5*y(i)^2/(2*a))./sqrt(a*2*pi);
        L2 = exp(-0.5*y(i)^2/(2*b))./sqrt(b*2*pi);
       if L1/L2 > tau
            H1 = H1+1;
        else
            H0 = H0+1;
        end
    end
    P_f = H1/N;
end