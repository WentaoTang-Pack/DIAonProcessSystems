clear;

%% sys
load A
load B
n = size(A,1);
measure_ind = [1,2,3,4];k_max = 2;
[C,Vn] = f_C(n, measure_ind);  %changing this matters 
Vd = eye(size(C,1));       % changing this matters
%% DARE to get K in u = Kx
Q = eye(n)*10; % Q and R in paper
R = eye(n)*1;
[P,~,~] = dare(A', C', Q, R); %% both K in two methods are u = -Kx
K = inv(R+B'*P*B)*B'*P*A;
L = P*C'/(C*P*C'+R);
K = - K;
%% Kf and covariance matrix
A_cl = A + B*K; %  stability
eig(A_cl);
[Sigma_xx, Sigma_xhatx,Sigma_hatxx,Sigma_hatxhatx,Sigma_yy,Sigma_xixi] = cov_matrix(A,B,K,L,C,Vd,Vn,zeros(n,n));
%% constraints on lambda
[M2] = f_M2(A,B,K,L,C,Sigma_xx,Sigma_hatxhatx,Sigma_xhatx,Sigma_hatxx);
Sigma_aa = zeros(n,n);
A_set = [];
lambda = 5;
for k = 1:k_max-1
    sensor_candi = setdiff(measure_ind,A_set);
    [ind] = f_bestSensor_ksparse(sensor_candi,Sigma_yy,M2,Sigma_aa,Sigma_xixi,A,B,K,L,C,Vd,Vn,lambda);  
    v(ind) = f_v_ksparse(Sigma_yy,M2,Sigma_aa,ind,lambda);
    temp = zeros(n,n); temp(ind,ind) = v(ind); % v is for different lambda
    Sigma_aa = Sigma_aa + temp;
    A_set = [A_set; ind];
end
sensor_candi = setdiff(measure_ind,A_set);
[ind] = f_bestSensor_ksparse(sensor_candi,Sigma_yy,M2,Sigma_aa,Sigma_xixi,A,B,K,L,C,Vd,Vn,lambda);  
[lambda_min,lambda_max] = f_lambda_lim_ksparse(Sigma_yy,M2,Sigma_aa,ind);
lambda_list = lambda_min+0.01:0.1:lambda_max - 0.01;
v_list = f_v_ksparse(Sigma_yy,M2,Sigma_aa,ind,lambda_list);

for i = 1:size(v_list,2)
    temp = zeros(n,n); temp(ind,ind) = v_list(i); % v is for different lambda
    Sigma_aa_temp = Sigma_aa + temp;
    [~,D1_k2(i),D2_k2(i)] = f_obj_ksparse(Sigma_xixi,Sigma_yy,Sigma_aa_temp,lambda,A,B,K,L,C,Vd,Vn);
    [P_d_k2(i),P_f_k2(i)] = f_AttackDetection(Sigma_yy,Sigma_aa_temp);
end   
plot(D2_k2,-D1_k2,'-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 1]), hold on;
save('D1_k2','D1_k2')
save('D2_k2','D2_k2')
save('P_d_k2','P_d_k2')
save('P_f_k2','P_f_k2')
% plot([D2(1),D2(i)],[-D1(1),D1(i)]) check if it is a straight line
% xlabel('$D(P_{y_a y_a} || P_{yy}$) ', 'Interpreter', 'latex', 'FontSize', 20);
% ylabel('$-D(P_{\xi_a \xi_a} || P_{\xi \xi}$) ', 'Interpreter', 'latex', 'FontSize', 20);
% grid on;
% legend('k = 1', 'k = 2', 'k =  3', 'k = 4', 'Location', 'northeast');
% xlim([-20 290])
% ylim([-630 85])
% annotation('arrow', [0.4, 0.6], [0.5, 0.8], 'LineWidth', 1.5); % Adjust coordinates
% 
% % Add LaTeX text (e.g., \lambda)
% text(pi, 0, '$\lambda_{\min}$', 'Interpreter', 'latex', 'FontSize', 14);