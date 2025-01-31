clear;
clc
close all
%% sys
load A
load B
n = size(A,1);
for j = 1:n
    measure_ind = [j];
    [C,Vn] = f_C(n, measure_ind);  %changing this matters 
    Vd = eye(size(C,1))*10*0.1;       % changing this matters
    %% DARE to get K in u = Kx
    Q = eye(n)*1; % Q and R in paper
    R = eye(n)*1;
    [P,~,~] = dare(A', C', Q, R); %% both K in two methods are u = -Kx
    K = inv(R+B'*P*B)*B'*P*A;
    L = P*C'/(C*P*C'+R);
    K = - K;
    %% Kf and covariance matrix
    A_cl = A + B*K; %stability
    eig(A_cl);
    [Sigma_xx, Sigma_xhatx,Sigma_hatxx,Sigma_hatxhatx,Sigma_yy,Sigma_xixi] = cov_matrix(A,B,K,L,C,Vd,Vn,zeros(n,n));
    %% constraints on lambda

    [M2] = f_M2(A,B,K,L,C,Sigma_xx,Sigma_hatxhatx,Sigma_xhatx,Sigma_hatxx,Sigma_xixi);
    [lambda_min,lambda_max] = f_lambda_lim_single(Sigma_yy,M2,measure_ind);
    for i = 1:lambda_max - lambda_min - 1 
        Sigma_aa = zeros(n,n);
        lambda = round(lambda_min) + i;
        lambda_list(i,j) = lambda;
        v(i,j) = f_v(Sigma_yy,M2,measure_ind,lambda); 
        Sigma_aa(measure_ind,measure_ind) = v(i,j);
        [~,~,~,Sigma_hatxahatxa,~,Sigma_xiaxia] = cov_matrix(A,B,K,L,C,Vd,Vn,Sigma_aa);
        [obj(i,j),D1(i,j),D2(i,j),D3(i,j)] = f_obj_s(Sigma_xixi,Sigma_xiaxia,Sigma_yy,Sigma_aa,lambda,measure_ind);
        [P_d(i,j),P_f(i,j)] = f_AttackDetection_s(Sigma_yy,v(i,j),measure_ind);
    end
end
 
%%
% plot(P_f,P_d), hold on;
yyaxis left;
plot(lambda_list(:,1), -D1(:,1)', '-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0 1]), hold on;
plot(lambda_list(1:19,2), -D1(1:19,2)', '-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0 1]), hold on;
plot(lambda_list(1:14,3), -D1(1:14,3)', '-^', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0 1]), hold on;
plot(lambda_list(1:14,4), -D1(1:14,4)', '-d', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0 1]), hold on;
ylim([-11 4])
xlim([4.3 24]);
ylabel('$-D(P_{\xi_a} || P_{\xi}$) ', 'Interpreter', 'latex', 'FontSize', 20);
ax = gca; % Get current axesfig
ax.YColor = [0 0 1];
ax.YLabel.FontName = 'Arial';
yyaxis right;
plot(lambda_list(:,1), P_d(:,1)', '-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
plot(lambda_list(1:19,2), P_d(1:19,2)', '-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
plot(lambda_list(1:14,3), P_d(1:14,3)', '-^', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
plot(lambda_list(1:14,4), P_d(1:14,4)', '-d', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
ylim([0.3 0.95]);
ylabel('Probability of detection', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
ax = gca; % Get current axes
ax.YColor = [1 0 0];
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Arial')
% title('Plot with Left and Right Y-Axes');
grid on;                   % Turn on the grid
legend('Attack on measurement 1', 'Attack on measurement 2', 'Attack on measurement 3', 'Attack on measurement 4','Attack on measurement 1', 'Attack on measurement 2', 'Attack on measurement 3', 'Attack on measurement 4', 'Location', 'northeast');

 



figure(2)
plot(D2(:,1),-D1(:,1),'-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [1 0 0]), hold on;
plot(D2(:,2),-D1(:,2),'-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 1]), hold on;
plot(D2(:,3),-D1(:,3),'-^', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [1 0.5 0]), hold on;
plot(D2(:,4),-D1(:,4),'-d', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.5 0.5 0.5]), hold on;
% plot([D2(1),D2(i)],[-D1(1),D1(i)])
% Adding labels and title
xlabel('$D(P_{y_a} || P_{y}$) ', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$-D(P_{\xi_a} || P_{\xi}$) ', 'Interpreter', 'latex', 'FontSize', 20);
% title('Mutual Information vs. KL Divergence');
% Adding grid and legend
grid on;
legend('Attack on measurement 1', 'Attack on measurement 2', 'Attack on measurement 3', 'Attack on measurement 4', 'Location', 'northeast');
annotation('arrow', [0.4, 0.6], [0.5, 0.8], 'LineWidth', 1.5); % Adjust coordinates
% 
% % Add LaTeX text (e.g., \lambda)
text(1, -2, '$\lambda_{\min}$', 'Interpreter', 'latex', 'FontSize', 14,'FontWeight', 'bold'); 

text(1, -2, '$\lambda_{\max}$', 'Interpreter', 'latex', 'FontSize', 14,'FontWeight', 'bold'); 


 

figure(3)
yyaxis left;
plot(lambda_list(:,1), D3(:,1)', '-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [0 0 1]), hold on;
ylim([-1 7])
xlim([4.3 25]);
ylabel('$D(P_{\hat{x}_a} || P_{\hat{x} }$) ', 'Interpreter', 'latex', 'FontSize', 20);
ax = gca; % Get current axes
ax.YColor = [0 0 1];
ax.YLabel.FontName = 'Arial';
yyaxis right;
plot(lambda_list(:,1),P_d(:,1), '-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
% plot(lambda_list(1:18,2), P_d(1:18,2)', '-s', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
% plot(lambda_list(1:15,3), P_d(1:15,3)', '-^', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
% plot(lambda_list(1:14,4), P_d(1:14,4)', '-d', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
ylim([0.3 0.95]);
ylabel('Probability of detection', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
ax = gca; % Get current axes
ax.YColor = [1 0 0];
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Arial')
% title('Plot with Left and Right Y-Axes');
grid on;                   % Turn on the grid
 
 