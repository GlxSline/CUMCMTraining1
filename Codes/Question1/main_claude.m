% r_0 = 8.8;
% r0_init = 8.8;
% k = 0.55 / (2 * pi);

% t = linspace(0, 300, 301);

% f = @(r) (r / 2) * ((r ^ 2 + k ^ 2) ^ 0.5 - (r_0 ^ 2 + k ^ 2) ^ 0.5) + ...
%     (0.5 * k * k) * ((log(r + (r ^ 2 + k ^ 2) ^ 0.5)) - log(r_0 + (r_0 ^ 2 + k ^ 2) ^ 0.5));

% r_sol = zeros(size(t));
% opts = optimset('Display', 'off');

% for i = 1:numel(t)
%     tk = t(i);
%     fun = @(rr) f(rr) + tk * k;
%     r_k = fzero(fun, r0_init, opts);
%     r_sol(i) = r_k;
%     r0_init = r_k;
% end

% theta = (r_sol - r_0) ./ k;

% vn = zeros(1, 301);
% sn = zeros(1, 301);
% for i = 1:numel(t)
%     vn(i) = ((k*k)/(k*k + r_sol(i)^2))^0.5;
%     if i>1
%         sn(i) = sn(i-1) + vn(i);
%     end
% end
% sn = 8.8-sn;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% r_0 = 8.8;
% r0_init = 8.8;
% k = 0.55 / (2 * pi);
% t = linspace(0, 300, 301);

% f1 = @(t) sqrt((8.8 + k.*t).^2 + k.^2);
% f = @(x) arrayfun(@(xi) integral(f1,0,xi), x);

% % rI = zeros(size(t));
% % sI = linspace(0, -32 * pi, 301);
% r_sol = zeros(size(t));
% opts = optimset('Display', 'off');

% for i = 1:numel(t)
%     ti = t(i);
%     fun = @(rr) f(rr) + ti * k;
%     r_k = fzero(fun, r0_init, opts);
%     r_sol(i) = r_k;
%     r0_init = r_k;
% end

% theta = (r_sol - r_0) ./ k;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% r0 = 8.8; % 初始半径
% k = 0.55 / (2 * pi); % 螺旋参数
% v = 1; % 运动速度
% t = linspace(0, 300, 301);

% % 局部弧长微分 d s / dθ
% dsdtheta = @(th) sqrt((r0 + k .* th) .^ 2 + k ^ 2);

% % 弧长函数 s(θ)=∫₀^θ dsdtheta(τ) dτ
% s_of_theta = @(th) arrayfun(@(x) integral(dsdtheta, 0, x), th);

% r_sol = zeros(size(t));
% theta0 = 0; % fzero 初猜，从 0 向负方向
% opts = optimset('Display', 'off');

% for i = 1:numel(t)
%     ti = t(i);
%     % 解 s(θ) + v*t = 0
%     fun = @(theta) s_of_theta(theta) + v * ti;
%     theta_i = fzero(fun, theta0, opts);
%     theta0 = theta_i; % 下次初猜用上次解
%     r_sol(i) = r0 + k * theta_i; % 计算 r(t)
% end

% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% fun = @(rho1, theta1, rho2, theta2, l) ...
%     [rho2 - rho1 - k * (theta2 - theta1); ...
%      rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2];

% numb = 20;
% result_rho = zeros(301, numb);
% result_theta = zeros(301, numb);
% result_rho(:, 1) = r_sol';
% result_theta(:, 1) = (r_sol' - 8.8) ./ k;

% for i = 1:numel(t)

%     for j = 1:(numb - 1)
%         rho1 = result_rho(i, j);
%         theta1 = result_theta(i, j);

%         if j == 1
%             fun_solve = @(x) fun(rho1, theta1, x(1), x(2), 2.86);
%         else
%             fun_solve = @(x) fun(rho1, theta1, x(1), x(2), 1.65);
%         end
%         x0 = [0; 0];

%         options = optimoptions('fsolve', 'Display', 'off');
%         sol = fsolve(fun_solve, x0, options);
%         result_rho(i, j+1) = sol(1);
%         result_theta(i, j+1) = sol(2);

%     end
% end
