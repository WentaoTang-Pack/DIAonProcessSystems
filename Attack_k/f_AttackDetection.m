function [P_d,P_f] = f_AttackDetection(Sigma_yy,Sigma_aa)
    tau = 1;
    N = 20000;
    Sigma_yaya = (Sigma_yy + Sigma_aa + (Sigma_yy + Sigma_aa)') / 2;
    ya = mvnrnd(zeros(1,size(Sigma_yy,1)), Sigma_yaya, N);
    y = mvnrnd(zeros(size(Sigma_yy,1),1), (Sigma_yy+Sigma_yy')/2, N);
    k = size(ya,2);
    % Test on probability of detection
    H1 = 0;
    H0 = 0;
    for i = 1:N  
        L1 = exp(-0.5*ya(i,:)*inv(Sigma_yaya)*ya(i,:)')./sqrt((2*pi)^k*det(Sigma_yaya));
        L2 = exp(-0.5*ya(i,:)*inv(Sigma_yy)*ya(i,:)')./sqrt((2*pi)^k*det(Sigma_yy));
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
        L1 = exp(-0.5*y(i,:)*inv(Sigma_yaya)*y(i,:)')./sqrt((2*pi)^k*det(Sigma_yaya));
        L2 = exp(-0.5*y(i,:)*inv(Sigma_yy)*y(i,:)')./sqrt((2*pi)^k*det(Sigma_yy));
        if L1/L2 > tau
            H1 = H1+1;
        else
            H0 = H0+1;
        end
    end
    P_f = H1/N;
end