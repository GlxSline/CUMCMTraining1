% %* 约束方程
% function F = segment_eq_1(x, rho1, theta1, k, l)
%     rho2 = x(1);
%     theta2 = x(2);
%     F = [
%          rho2 - rho1 - k * (theta2 - theta1);
%          rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
%          ];

% end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% %* -100~0
% % r = kθ
% %* rho1
% width = 1.7;
% k = width / (2 * pi);
% v = 1;
% t = linspace(0, -100, 101);
% bench_numb = 224;
% r0 = 4.5;

% t_of_theta1 = @(theta) (theta * sqrt(theta ^ 2 + 1) / 2 + log(theta + sqrt(theta ^ 2 + 1)) / 2) * k;

% r_sol = zeros(numel(t), 1);
% theta_sol = zeros(numel(t), 1);
% opts_fz = optimset('Display', 'off');
% theta_start = r0 / k;
% t_offset = t_of_theta1(theta_start);

% for i = 1:numel(t)
%     ti = t(i);
%     funs = @(th) t_of_theta1(th) - v * ti - v * t_offset;
%     si = fzero(funs, theta0, opts_fz);
%     theta_sol(i) = si;
%     r_sol(i) = k * si;
%     theta0 = si;
% end

% result_rho = zeros(numel(t), bench_numb);
% result_theta = zeros(numel(t), bench_numb);
% result_rho(:, 1) = r_sol;
% result_theta(:, 1) = theta_sol;

%% ————————————————
%  参数设定
%% ————————————————
width      = 1.7;
k          = width/(2*pi);
v          = 1;
T          = 100;           % 总时间
N          = 101;           % 时间步数
t          = linspace(0, T, N)';  % 列向量
bench_numb = 224;

% 弧长→时间映射（原点起算）
t_of_theta = @(th) (k/(2*v)) * ( ...
    th.*sqrt(th.^2+1) + ...
    log(th + sqrt(th.^2+1)) ...
);

%% ————————————————
%  每圈初始半径 r0_j 逐渐增大
%  例：从 4.5 开始，每圈增加 0.1
%% ————————————————
r0_start = 4.5;
delta_r  = 0.1;
r0s      = r0_start + (0:bench_numb-1)' * delta_r;  % [224×1] 向量

% 对应的初始极角 θ0_j 和时间偏移 t_offset_j
theta0s    = r0s / k;                      % [224×1]
t_offsets  = t_of_theta(theta0s);          % [224×1]

%% ————————————————
%  数值求解 θ(i,j)、ρ(i,j)
%% ————————————————
theta_sol = zeros(N, bench_numb);
rho_sol   = zeros(N, bench_numb);

opts = optimset('Display','off');

for j = 1:bench_numb
    th_guess = theta0s(j);   % fzero 初值设为对应圈的起始角
    for i = 1:N
        % 方程：t_of_theta(th) - t_offsets(j) = t(i)
        F   = @(th) t_of_theta(th) - t_offsets(j) - t(i);
        thi = fzero(F, th_guess, opts);
        theta_sol(i,j) = thi;
        rho_sol(i,j)   = k * thi;
        th_guess       = thi; % 保持初值连续
    end
end

% 结果保存在 rho_sol、theta_sol 中


% %* rhoi

% options = optimoptions('fsolve', 'Display', 'off');

% for i = 1:numel(t)

%     for j = 1:(bench_numb - 1)
%         rho1 = result_rho(i, j);
%         theta1 = result_theta(i, j);
%         l = (j == 1) * 2.86 + (j > 1) * 1.65;
%         dth0 = l / rho1;
%         x0 = [rho1 + k * dth0; theta1 + dth0];
%         sol = fsolve(@(x) segment_eq_1(x, rho1, theta1, k, l), x0, options);
%         result_rho(i, j + 1) = sol(1);
%         result_theta(i, j + 1) = sol(2);
%     end

% end
