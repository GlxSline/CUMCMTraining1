theta1 = linspace(0, 32 * pi, 1000);
r1 = theta1 .* 1.7 ./ 2 ./ pi;
r2 = (theta1 + pi) .* 1.7 ./ 2 ./ pi;

x1 = r1 .* cos(theta1);
y1 = r1 .* sin(theta1);

x2 = r2 .* cos(theta1);
y2 = r2 .* sin(theta1);

theta = linspace(0, 2 * pi, 1000);
r = 4.5;
x3 = r .* cos(theta);
y3 = r .* sin(theta);

figure;
plot(x1, y1, 'LineWidth', 1.5);
axis equal; % 保持纵横比例相等
grid on; % 如需网格可打开，不需要可注释此行
hold on;

plot(x2, y2, 'LineWidth', 1.5);
axis equal; % 保持纵横比例相等
grid on; % 如需网格可打开，不需要可注释此行

plot(x3, y3, 'LineWidth', 1.5);
axis equal; % 保持纵横比例相等
grid on; % 如需网格可打开，不需要可注释此行
