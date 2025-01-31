% pick up the index for each measurement
% lambda = 8; -> 46th, 39th, 31th, 17th
D(1) = D1_k1(46);
D(2) = D1_k2(39);
D(3) = D1_k3(31);
D(4) = D1_k4(17);

P_d(1) = P_d_k1(46);
P_d(2) = P_d_k2(39);
P_d(3) = P_d_k3(31);
P_d(4) = P_d_k4(17);


P_f(1) = P_f_k1(46);
P_f(2) = P_f_k2(39);
P_f(3) = P_f_k3(31);
P_f(4) = P_f_k4(17);


yyaxis left;
plot(-D, '-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', 'b'), hold on;
ylim([-60 10])
xlim([0.8 4.2]);
ylabel('$-D(P_{\xi_a} || P_{\xi}$) ', 'Interpreter', 'latex', 'FontSize', 20);
ax = gca; % Get current axes
ax.YColor = [0 0 1];
ax.YLabel.FontName = 'Arial';
yyaxis right;
plot(P_d, '-o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
plot(P_f, '--o', 'LineWidth', 1.5, 'MarkerSize', 8, 'Color', [1 0 0]), hold on;
ylim([0 1.3]);
ylabel('Detection', 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
ax = gca; % Get current axes
ax.YColor = [1 0 0];
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 18, 'FontName', 'Arial')
% title('Plot with Left and Right Y-Axes');
grid on;                   % Turn on the grid
legend('$-D(P_{\xi_a} || P_{\xi}$)','Probability of detection','Probability of false alarm','Interpreter', 'latex');

 