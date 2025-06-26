% k = 0.55 / (2 * pi);
% r_0 = 8.8;

% f = @(r) (r / 2) * ((r ^ 2 + k ^ 2) ^ 0.5 - (r_0 ^ 2 + k ^ 2) ^ 0.5) + ...
%     (0.5 * k * k) * ((log(r + (r ^ 2 + k ^ 2) ^ 0.5)) - log(r_0 + (r_0 ^ 2 + k ^ 2) ^ 0.5));

% rr = linspace(0, 8.8, 1000);
% t = zeros(1, 1000);

% for i = 1:1000
%     t(i) = f(rr(i)) / (-0.55);
% end

% % 绘图
% figure;
% plot(t, rr, 'LineWidth', 1.5);
% grid on;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% f = @(th) (k / 2) * (th * (th * th + 1) ^ (1/2) + log(th + (th * th + 1) ^ (1/2)));
% th0 = 0;
% th1 = linspace(0, 32* pi, 1000);
% l = zeros(1, 1000);
% for i = 1:1000
%     l(i) = f(th1(i))- f(th0);
% end

% % 绘图
% figure;
% plot(th1, l, 'LineWidth', 1.5);
% grid on;

% % 标注
% title('等距螺旋线 (Archimedean Spiral)');
% xlabel('X 轴');
% ylabel('Y 轴');

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% k = 0.55/(2*pi);
% r_0 = 8.8;

% syms r
% f1 = (r^2 + k^2)^(1/2);
% F = int(f1, x);

% rr = linspace(0, 8.8, 1000);
% t = zeros(1, 1000);

% for i = 1:1000
%     t(i) = f(rr(i)) / (-0.55);
% end

% % 绘图
% figure;
% plot(t, rr, 'LineWidth', 1.5);
% grid on;

% % 标注
% title('等距螺旋线 (Archimedean Spiral)');
% xlabel('X 轴');
% ylabel('Y 轴');

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% k = 0.55;
% fr = @(rr) (rr / 2) * ((rr ^ 2 + k ^ 2) ^ (1/2)) + (k * k / 2) * log(rr + (rr ^ 2 + k ^ 2) ^ (1/2));

% r0 = 8.8;
% r1 = linspace(8.8, 0, 1000);
% t = zeros(1, 1000);

% for i = 1:1000
%     t(i) = (fr(r1(i)) - fr(r0)) / (-0.55);
% end

% figure;
% plot(r1, t, 'LineWidth', 1.5);
% grid on;

% xlabel('X 轴');
% ylabel('Y 轴');

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% k = 0.55 / (2 * pi);

% f = @(t) ((8.8 + k * t) .^ 2 + k * k) .^ (1/2);

% F = @(x) arrayfun(@(xx) integral(f, 0, xx), x);

% s = linspace(0, -32 * pi, 1000);
% t = zeros(1, 1000);
% r = zeros(1, 1000);

% for i = 1:1000
%     t(i) = F(s(i));
%     r(i) = 8.8 + k * s(i);
% end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% 参数
k = 0.55/(2*pi);

% 定义原函数
f = @(t) sqrt((8.8 + k*t).^2 + k^2);
F = @(x) arrayfun(@(xx) integral(f,0,xx), x);

% 在 s∈[s_min,s_max] 上预计算
s_vals = linspace(-32*pi, 0, 1000);    % 注意这里按从小到大排列
t_vals = F(s_vals);

% 构造反演函数 s = s_of_t(t)
s_of_t = @(tq) interp1(t_vals, s_vals, tq, 'pchip');

% 验证：把原来算得的 t 作为输入，重算 s
s = linspace(0, -32*pi, 1000);
t = arrayfun(@(ss) F(ss), s);
t_rec = linspace(0, 300, 301);
s_rec = s_of_t(t_rec);   % 应该近似等于原 s
r_rec = 8.8 - k .* s_rec;
x_rec = zeros(1, 301);
y_rec = zeros(1, 301);
r = 8.8 + k .* s;
x = zeros(1, 1000);
y = zeros(1, 1000);
s = s.*(-1);
t = t.*(-1);
for i = 1:301
    x_rec(i) = r_rec(i) * cos(s_rec(i));
    y_rec(i) = r_rec(i) * sin(s_rec(i));

end
for i = 1:1000
    x(i) = r(i) * cos(s(i));
    y(i) = r(i) * sin(s(i));

end



% 作图对比
figure;
plot(t, s, 'r', 'LineWidth', 0.5);
hold on;
plot(t_rec, s_rec,'b', 'LineWidth', 0.5);
xlabel('t');
ylabel('s');

legend('原 s(t)','反演 s_{rec}(t)');
grid on;
