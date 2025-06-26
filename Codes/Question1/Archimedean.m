a = 8.8;        % 起始半径
d = 0.55;        % 相邻两圈径向间距
b = d/(2*pi); % 螺距参数

theta = linspace(0, -32*pi, 2000); % θ 从 0 到 6π，共三圈
r     = a + b*theta;             % 极坐标半径

% 极坐标转直角坐标
x = r .* cos(theta);
y = r .* sin(theta);

% 绘图
figure;
plot(x, y, 'LineWidth', 1.5);
axis equal;       % 保持纵横比例相等
grid on;       % 如需网格可打开，不需要可注释此行

ax = gca;        % 获取当前坐标轴句柄
ax.XTick = [];   % 隐藏 X 轴刻度
ax.YTick = [];   % 隐藏 Y 轴刻度

