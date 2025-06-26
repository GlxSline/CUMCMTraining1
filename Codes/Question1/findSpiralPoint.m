function [x1, y1, theta1] = findSpiralPoint(x0, y0, d, a, b)
% findSpiralPoint  在 Archimedean 螺旋线上，已知一点 (x0,y0) 与
% 直线距离 d，求另一点 (x1,y1)。
%
% 输入：
%   x0, y0   — 已知点坐标
%   d        — 两点之间的直线距离
%   a, b     — 螺旋线参数：r = a + b*theta
% 输出：
%   x1, y1   — 求得的另一点坐标
%   theta1   — 求得点的极角

    % 1) 计算已知点的极角和极径
    theta0 = atan2(y0, x0);
    r0     = hypot(x0, y0);

    % 2) 构造待解方程 f(theta1)=0
    f = @(theta1) ...
        sqrt( ((a + b*theta1).*cos(theta1) - x0).^2 + ...
              ((a + b*theta1).*sin(theta1) - y0).^2 ) ...
        - d;

    % 3) 用 fzero 求根，初始猜测可取 theta0 + d/r0
    theta1 = fzero(f, theta0 + d/r0);

    % 4) 计算解点坐标
    r1 = a + b*theta1;
    x1 = r1 * cos(theta1);
    y1 = r1 * sin(theta1);
end
