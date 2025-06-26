%* 约束方程
function F = segment_eq_1(x, rho1, theta1, k, l)
    rho2 = x(1);
    theta2 = x(2);
    F = [
         rho2 - rho1 - k * (theta2 - theta1);
         rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
         ];

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% * -100~0
% r = kθ
%* rho1
width = 1.7;
k = width / (2 * pi);
v = 1;
t = linspace(0, 100, 101);
bench_numb = 224;
r0 = 4.5;

t_of_theta1 = @(theta) (theta * sqrt(theta ^ 2 + 1) / 2 + ...
    log(theta + sqrt(theta ^ 2 + 1)) / 2) * k;

r_sol = zeros(numel(t), 1);
theta_sol = zeros(numel(t), 1);
opts_fz = optimset('Display', 'off');
theta0 = r0 / k;
t_offset = t_of_theta1(theta0);

for i = 1:numel(t)
    ti = t(i);
    funs = @(th) t_of_theta1(th) - v * ti - v * t_offset;
    si = fzero(funs, theta0, opts_fz);
    theta_sol(i) = si;
    r_sol(i) = k * si;
    theta0 = si + 0.2;
end

%* rhoi

result_rho = zeros(numel(t), bench_numb);
result_theta = zeros(numel(t), bench_numb);
result_rho(:, 1) = r_sol;
result_theta(:, 1) = theta_sol;

options = optimoptions('fsolve', 'Display', 'off');

for i = 1:numel(t)

    for j = 1:(bench_numb - 1)
        rho1 = result_rho(i, j);
        theta1 = result_theta(i, j);
        l = (j == 1) * 2.86 + (j > 1) * 1.65;
        % dth0 = l / rho1;
        % x0 = [rho1 + k * dth0; theta1 + dth0];
        x0 = [result_rho(i, j), result_theta(i, j)];
        sol = fsolve(@(x) segment_eq_1(x, rho1, theta1, k, l), x0, options);
        result_rho(i, j + 1) = sol(1);
        result_theta(i, j + 1) = sol(2);
    end

end

result_ki = zeros(numel(t), bench_numb);
result_x = zeros(numel(t), bench_numb);
result_y = zeros(numel(t), bench_numb);
result_alpha = zeros(numel(t), bench_numb - 1);
result_beta = zeros(numel(t), bench_numb - 1);
result_k = zeros(numel(t), bench_numb - 1);
result_v = zeros(numel(t), bench_numb);

for i = 1:numel(t)

    for j = 1:bench_numb
        th = result_theta(i, j);
        result_ki(i, j) = ((8.8 + k * th) * cos(th) + k * sin(th)) / (k * cos(th) - (k * th + 8.8) * sin(th));
        result_x(i, j) = result_rho(i, j) * cos(th);
        result_y(i, j) = result_rho(i, j) * sin(th);

    end

end

for i = 1:numel(t)

    for j = 1:(bench_numb - 1)
        result_k(i, j) = (result_y(i, j + 1) - result_y(i, j)) / (result_x(i, j + 1) - result_x(i, j));
    end

end

for i = 1:numel(t)

    for j = 1:(bench_numb - 1)
        result_alpha(i, j) = atan(abs((result_ki(i, j) - result_k(i, j)) / (result_ki(i, j) * result_k(i, j) +1)));
        result_beta(i, j) = atan(abs((result_ki(i, j + 1) - result_k(i, j)) / (result_ki(i, j + 1) * result_k(i, j) +1)));
    end

    result_v(i, 1) = 1;

    for j = 1:(bench_numb - 1)
        result_v(i, j + 1) = result_v(i, j) * cos(result_alpha(i, j)) / cos(result_beta(i, j));
    end

end