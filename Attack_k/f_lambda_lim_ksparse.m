function [lambda_min,lambda_max] = f_lambda_lim_ksparse(Sigma_yy,M2,Sigma_aa,i)
        n = size(Sigma_yy,2);
        C = inv(Sigma_yy + Sigma_aa);
        D = inv(eye(n)+M2*inv(Sigma_yy))*M2;
        %% solve inequalities for lambda
        syms x
        a = 1/Sigma_yy(i,i);
        b = M2(i,i);
        c = C(i,i);
        d = D(i,i);
        inequality1 = a * x - b > 0;
        inequality2 = ((a*x-b)*(c+d) +c*d*(1-x))^2 -4*(a*x-b)*c*d*(d-c*x) > 0; % v is real
        inequality3 = (b - a*x)*(c+d) + c*d*(x-1) > 0; % v is positive
        solution = solve([inequality1,inequality2, inequality3], x, 'ReturnConditions', true);  % Solve the inequality
        cond = solution.conditions;
        cell = num2cell(vpa(cond,10));
        ineq = cell{1};
        bounds = children(ineq); 
        if lhs(bounds{1}) == x
            lambda_max = double(rhs(bounds{1}));  
            lambda_min = double(lhs(bounds{2})); 
        else
            lambda_min = double(lhs(bounds{1}));  
            lambda_max = double(rhs(bounds{2}));  

        end
        disp(lambda_min);
        disp(lambda_max);
end