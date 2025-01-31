clear;
clc
close all
%% sys
load A
load B
n = size(A,1);
j = 2; % choose to attack the first measurement as an example
    measure_ind = j;
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
    i = 3;
    Sigma_aa = zeros(n,n);
        lambda = round(lambda_min) + i;
        lambda_list(i) = lambda; % use the smallest lambda to make the Gaussian more obvious
        v(i,j) = f_v(Sigma_yy,M2,measure_ind,lambda); 
        Sigma_aa(measure_ind,measure_ind) = v(i,j);
    [Sigma_xaxa, ~,~,Sigma_hatxahatxa,~,~] = cov_matrix(A,B,K,L,C,Vd,Vn,Sigma_aa);
    
  
 