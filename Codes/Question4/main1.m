% %* 龙头坐标方程
% function [rho, theta, x, y] = t_of_theta2(ta)
%     % x = rho * cos (theta);
%     % y = rho * sin (theta);
%     % r_AE1 = sqrt((x - x_E1) ^ 2 + (y - y_E1) ^ 2);
%     % k_AE1 = (y - y_E1) / (x - x_E1);

%     if ta < t1 %* IV
%         [rho, theta, x, y] = calculate2(ta);
%     elseif ta < t2 %* II
%         [rho, theta, x, y] = calculate3(ta);
%     else %* III

%     end

% end

% function [rho, theta, x, y] = calculate2(ta)
%     theta_AE1 = v * ta / r_E1E2;
%     vector_E2E1 = [x_E1 - x_E2; y_E1 - y_E2];
%     rotation_matrix = [cos(theta_AE1), sin(theta_AE1); -1 * sin(theta_AE1), cos(theta_AE1)];
%     vector_E2A = rotation_matrix * vector_E2E1;
%     x = vector_E2A(1) + x_E2;
%     y = vector_E2A(2) + y_E2;
%     rho = sqrt(x ^ 2 + y ^ 2);

%     if x > 0

%         if y >= 0
%             theta = asin(y / rho);

%         else
%             theta = pi - asin(y / rho);
%         end

%     elseif x <= 0

%         if y >= 0
%             theta = asin(y / rho);
%         else
%             theta = -1 * pi - asin (y / rho);
%         end

%     end

% end

% function [rho, theta, x, y] = calculate3(ta)
%     theta_AE3 = v * (ta - t1) / r_E3E4;
%     vector_E4E3 = [x_E3 - x_E4; y_E3 - y_E4];
%     rotation_matrix = [cos(theta_AE3), -1 * sin(theta_AE3); sin(theta_AE3), cos(theta_AE3)];
%     vector_E4A = rotation_matrix * vector_E4E3;
%     x = vector_E4A(1) + x_E4;
%     y = vector_E4A(2) + y_E4;
%     rho = sqrt(x ^ 2 + y ^ 2);

%     if x > 0

%         if y >= 0
%             theta = asin(y / rho);

%         else
%             theta = pi - asin(y / rho);
%         end

%     elseif x <= 0

%         if y >= 0
%             theta = asin(y / rho);
%         else
%             theta = -1 * pi - asin (y / rho);
%         end

%     end
% end

% function [rho, theta, x, y] = calculate4(ta)
%     %* r = k (θ + Π/2)；
    

% end

% function F = theta2_of_theta1(input)

% end

% %* 约束方程
% function F = segment_eq_1(x, rho1, theta1, k, l)
%     rho2 = x(1);
%     theta2 = x(2);
%     F = [
%          rho2 - rho1 - k * (theta2 - theta1);
%          rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
%          ];

% end

% function F = segment_eq_2(x, rho1, theta1, R, xE, yE, l)
%     rho2 = x(1);
%     theta2 = x(2);
%     x1 = rho1 * cos (theta1);
%     x2 = rho2 * cos (theta2);
%     y1 = rho1 * sin (theta1);
%     y2 = rho2 * sin (theta2);
%     k1 = (y1 - yE) / (x1 - xE);
%     k2 = (y2 - yE) / (x2 - xE);

%     vector1 = [x1 - xE; y1 - yE];
%     vector2 = [x2 - xE; y2 - yE];
%     % l1 = sqrt((y1 - yE)^2 + (x1 - xE)^2);
%     % l2 = sqrt((y2 - yE)^2 + (x2 - xE)^2);
%     theta = atan((k1 - k2) / (1 + k1 * k2));

%     rotation_matrix = [cos(theta), sin(theta); -1 * sin(theta), cos(theta)];

%     F = [
%          vector2 - rotation_matrix * vector1;
%          R ^ 2 + R ^ 2 - 2 * R * R * cos(theta) - l ^ 2;
%          ];

% end

% %* 计算x_E, y_E
% theta_E = 4.5 * 2 * pi / 1.7;
% rho_E = 4.5;
% x_E1 = rho_E * cos(theta_E);
% y_E1 = rho_E * sin(theta_E);
% x_E3 = x_E1 / (-3);
% y_E3 = y_E1 / (-3);
% x_E5 = -1 * x_E1;
% y_E5 = -1 * y_E1;

% k_E1 = (sin(theta_E) + theta_E * cos(theta_E)) / (cos(theta_E) - theta_E * sin(theta_E));
% k_E1E5 = y_E1 / x_E1;
% k_E1E2 = -1 / k_E1;
% theta_E2 = atan(k_E1E2);

% gamma_E = atan((k_E1E5 - k_E1E2) / (1 + k_E1E2 * k_E1E5));
% r_E1E2 = 3 / cos(gamma_E);
% r_E3E4 = r_E1E2 / 2;
% r_E1E3 = sqrt((x_E1 - x_E3) ^ 2 + (y_E1 - y_E3) ^ 2);

% x_E2 = x_E1 + r_E1E2 * cos(theta_E2);
% y_E2 = y_E1 + r_E1E2 * sin(theta_E2);
% x_E4 = (3 * x_E3 - x_E2) / 2;
% y_E4 = (3 * y_E3 - y_E2) / 2;

% k_E2E3 = (y3 - y2) / (x3 - x2);

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%* -100~0
% r = kθ
%* rho1
width = 1.7;
k = width / (2 * pi);
v = 1;
t = linspace(0, -100, 101);
bench_numb = 224;

t_of_theta1 = @(theta) theta * sqrt(theta ^ 2 + 1) / 2 + log(theta + sqrt(theta ^ 2 + 1));

r_sol = zeros(size(t));
theta_sol = zeros(size(t));
opts_fz = optimset('Display', 'off');
theta0 = 0;

for i = 1:numel(t)
    ti = t(i);
    funs = @(th) t_of_theta1(th) - v * ti;
    si = fzero(funs, theta0, opts_fz);
    theta_sol(i) = si;
    r_sol(i) = k * si;
    theta0 = si;
end

result_rho = zeros(numel(t), numb);
result_theta = zeros(numel(t), numb);
result_rho(:, 1) = r_sol;
result_theta(:, 1) = theta_sol;

%* rhoi

options = optimoptions('fsolve', 'Display', 'off');

for i = 1:numel(t)

    for j = 1:(numb - 1)
        rho1 = result_rho(i, j);
        theta1 = result_theta(i, j);
        l = (j == 1) * 2.86 + (j > 1) * 1.65;
        dth0 = l / rho1;
        x0 = [rho1 + k * dth0; theta1 + dth0];
        sol = fsolve(@(x) segment_eq_1(x, rho1, theta1, k, l), x0, options);
        result_rho(i, j + 1) = sol(1);
        result_theta(i, j + 1) = sol(2);
    end

end

% result_ki = zeros(numel(t), numb);
% result_x = zeros(numel(t), numb);
% result_y = zeros(numel(t), numb);
% result_alpha = zeros(numel(t), numb - 1);
% result_beta = zeros(numel(t), numb - 1);
% result_k = zeros(numel(t), numb - 1);
% result_v = zeros(numel(t), numb);

% for i = 1:numel(t)

%     for j = 1:numb
%         th = result_theta(i, j);
%         result_ki(i, j) = ((8.8 + k * th) * cos(th) + k * sin(th)) / (k * cos(th) - (k * th + 8.8) * sin(th));
%         result_x(i, j) = result_rho(i, j) * cos(th);
%         result_y(i, j) = result_rho(i, j) * sin(th);

%     end

% end

% for i = 1:numel(t)

%     for j = 1:(numb - 1)
%         result_k(i, j) = (result_y(i, j + 1) - result_y(i, j)) / (result_x(i, j + 1) - result_x(i, j));
%     end

% end

% for i = 1:numel(t)

%     for j = 1:(numb - 1)
%         result_alpha(i, j) = atan(abs((result_ki(i, j) - result_k(i, j)) / (result_ki(i, j) * result_k(i, j) +1)));
%         result_beta(i, j) = atan(abs((result_ki(i, j + 1) - result_k(i, j)) / (result_ki(i, j + 1) * result_k(i, j) +1)));
%     end

%     result_v(i, 1) = 1;

%     for j = 1:(numb - 1)
%         result_v(i, j + 1) = result_v(i, j) * cos(result_alpha(i, j)) / cos(result_beta(i, j));
%     end

% end

% %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% %* 0~100

% t = linspace(0, 100, 101);

% %* t1
% theta_E1E3 = atan((k_E1E2 - k_E2E3) / (1 + k_E1E2 * k_E2E3));
% l_E1E3 = theta_E1E3 * r_E1E2;
% t1 = l_E1E3 / v;

% %* t2
% l_E3E5 = theta_E1E3 * r_E3E4;
% t2 = t1 + l_E3E5 / v;

% %* x1, y1, rho1, theta1
% result_rho