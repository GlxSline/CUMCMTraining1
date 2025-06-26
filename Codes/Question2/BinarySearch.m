t_first = 400;
t_last = 420;

% for binary_times = 1:30
r0 = 8.8;
k = 0.55 / (2 * pi);
v = 1;
numb = 224;

t_numb = 41;

% if binary_times ~= 1
% t_first = t_last - (t_last - t_first) / t_numb;
% end

t = linspace(t_first, t_last, t_numb).';
disp(t_first);
disp(t_last);

A = @(u) u .* sqrt(u .^ 2 + k ^ 2) + k ^ 2 * log(u + sqrt(u .^ 2 + k ^ 2));
s_of_theta = @(theta) (A(r0 + k * theta) - A(r0)) / (2 * k);

r_sol = zeros(size(t));
theta_sol = zeros(size(t));
opts_fz = optimset('Display', 'off');
theta0 = 0;

for i = 1:t_numb
    ti = t(i);
    funs = @(th) s_of_theta(th) + v * ti;
    si = fzero(funs, theta0, opts_fz);
    theta_sol(i) = si;
    r_sol(i) = r0 + k * si;
    theta0 = si;
end

result_rho = zeros(t_numb, numb);
result_theta = zeros(t_numb, numb);
result_rho(:, 1) = r_sol;
result_theta(:, 1) = theta_sol;

options = optimoptions('fsolve', 'Display', 'off');

for i = 1:t_numb

    for j = 1:(numb - 1)
        rho1 = result_rho(i, j);
        theta1 = result_theta(i, j);
        l = (j == 1) * 2.86 + (j > 1) * 1.65;
        dth0 = l / rho1;
        x0 = [rho1 + k * dth0; theta1 + dth0];
        sol = fsolve(@(x) segment_eq(x, rho1, theta1, k, l), x0, options);
        result_rho(i, j + 1) = sol(1);
        result_theta(i, j + 1) = sol(2);
    end

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

result_ki = zeros(t_numb, numb);
result_x = zeros(t_numb, numb);
result_y = zeros(t_numb, numb);
result_alpha = zeros(t_numb, numb - 1);
result_beta = zeros(t_numb, numb - 1);
result_k = zeros(t_numb, numb - 1);
result_v = zeros(t_numb, numb);

for i = 1:t_numb

    for j = 1:numb
        th = result_theta(i, j);
        result_ki(i, j) = ((8.8 + k * th) * cos(th) + k * sin(th)) / (k * cos(th) - (k * th + 8.8) * sin(th));
        result_x(i, j) = result_rho(i, j) * cos(th);
        result_y(i, j) = result_rho(i, j) * sin(th);

    end

end

for i = 1:t_numb

    for j = 1:(numb - 1)
        result_k(i, j) = (result_y(i, j + 1) - result_y(i, j)) / (result_x(i, j + 1) - result_x(i, j));
    end

end

for i = 1:t_numb

    for j = 1:(numb - 1)
        result_alpha(i, j) = atan(abs((result_ki(i, j) - result_k(i, j)) / (result_ki(i, j) * result_k(i, j) +1)));
        result_beta(i, j) = atan(abs((result_ki(i, j + 1) - result_k(i, j)) / (result_ki(i, j + 1) * result_k(i, j) +1)));
    end

    result_v(i, 1) = 1;

    for j = 1:(numb - 1)
        result_v(i, j + 1) = result_v(i, j) * cos(result_alpha(i, j)) / cos(result_beta(i, j));
    end

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
result_b = zeros(t_numb, numb - 1);
result_k1 = zeros(t_numb, 1);
result_k2 = zeros(t_numb, 1);
result_x1 = zeros(t_numb, 1);
result_y1 = zeros(t_numb, 1);
result_x2 = zeros(t_numb, 1);
result_y2 = zeros(t_numb, 1);

for i = 1:t_numb

    for j = 1:numb - 1
        result_b(i, j) = result_y(i, j) - result_k(i, j) * result_x(i, j);

    end

    result_k1(i, 1) = (0.275 * result_k(i, 1) + 0.15) / (0.275 - 0.15 * result_k(i, 1));
    result_k2(i, 1) = (0.275 * result_k(i, 1) - 0.15) / (0.275 + 0.15 * result_k(i, 1));

end

% 计算x1, y1

fun_x1_y1 = @(x1, y1, k, b, k1, xb, yb) ...
    [ ...
     y1 - k1 * (x1 - xb) - yb; ...
     abs(k * x1 + b - y1) - 0.15 * sqrt(k ^ 2 + 1)];

options = optimoptions('fsolve', 'Display', 'off');
init_guess = [0, 0];
x1_y1_sol = zeros(1, 2);

for i = 1:t_numb
    k = result_k(i, 1);
    k1 = result_k1(i, 1);
    b = result_b(i, 1);
    xb = result_x(i, 1);
    yb = result_y(i, 1);
    fun_solve = @(xy) fun_x1_y1(xy(1), xy(2), k, b, k1, xb, yb);

    x1_y1 = fsolve(fun_solve, init_guess, options);
    x1_y1_sol(1) = x1_y1(1); % x1
    x1_y1_sol(2) = x1_y1(2); % y1

    init_guess = x1_y1_sol;
    result_x1(i, 1) = x1_y1_sol(1);
    result_y1(i, 1) = x1_y1_sol(2);

end

% 计算x2, y2

fun_x2_y2 = @(x2, y2, k, b, k2, xb, yb) ...
    [ ...
     y2 - k2 * (x2 - xb) - yb; ...
     (k * x2 + b - y2) ^ 2 - (0.15 * sqrt(k ^ 2 + 1)) ^ 2];

options = optimoptions('fsolve', 'Display', 'off');
init_guess = [0, 0];
x2_y2_sol = zeros(1, 2);

for i = 1:t_numb
    k = result_k(i, 1);
    k2 = result_k2(i, 1);
    b = result_b(i, 1);
    xb = result_x(i, 1);
    yb = result_y(i, 1);
    fun_solve = @(xy) fun_x2_y2(xy(1), xy(2), k, b, k2, xb, yb);

    x2_y2 = fsolve(fun_solve, init_guess, options);
    x2_y2_sol(1) = x2_y2(1); % x2
    x2_y2_sol(2) = x2_y2(2); % y2

    init_guess = x2_y2_sol;
    result_x2(i, 1) = x2_y2_sol(1);
    result_y2(i, 1) = x2_y2_sol(2);

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% 计算每个x1, y1到板凳的距离

result_d = zeros(t_numb, numb - 3, 2);

for i = 1:t_numb

    for j = 1:numb - 3
        x = [result_x1(i, 1), result_x2(i, 1)];
        y = [result_y1(i, 1), result_y2(i, 1)];

        for ii = 1:2
            result_d(i, j, ii) = abs((result_k(i, j + 2) * x(ii) - y(ii) + result_b(i, j + 2))) / sqrt(result_k(i, j + 2) ^ 2 +1);
        end

    end

end

min_d1 = min(result_d(:, :, 1), [], 2);
min_d2 = min(result_d(:, :, 2), [], 2);
min_d = [min_d1, min_d2];

for i = 1:t_numb

    if min_d1(i) < 0.15
        break;
    elseif min_d2(i) < 0.15
        break;
    end

end

disp(i);
t_last = t_first + (t_last - t_first) / (((t_numb - 1))) * (i - 1);
disp(t_last);

% end

function F = segment_eq(x, rho1, theta1, k, l)
    rho2 = x(1);
    theta2 = x(2);
    F = [
         rho2 - rho1 - k * (theta2 - theta1);
         rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
         ];

end
