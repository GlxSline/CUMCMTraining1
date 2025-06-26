t_first = 411;
t_last = 413;
r0 = 8.8;
k = 0.55 / (2 * pi);
v = 1;
bench_numb = 224;

result_d_all = zeros(bench_numb - 3, 2, 30);

for binary_times = 1:30
    t_mid = (t_first + t_last) / 2;
    A = @(u) u .* sqrt(u .^ 2 + k ^ 2) + k ^ 2 * log(u + sqrt(u .^ 2 + k ^ 2));
    s_of_theta = @(theta) (A(r0 + k * theta) - A(r0)) / (2 * k) + v * t_mid;
    opts_fz = optimset( ...
        'Display', 'off', ...
        'TolX', 1e-10, ...
        'TolFun', 1e-10, ...
        'MaxIter', 1000, ...
        'MaxFunEvals', 1000);
    theta0 = -1 * v * t_mid / r0;
    theta_of_t = fzero(s_of_theta, theta0, opts_fz);
    rho_of_t = r0 + k * theta_of_t;

    bench_rho = zeros(1, bench_numb);
    bench_theta = zeros(1, bench_numb);
    bench_rho(1) = rho_of_t;
    bench_theta(1) = theta_of_t;

    options = optimoptions( ...
        'fsolve', ...
        'Display', 'off', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 1000, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);

    for i = 1:bench_numb - 1
        rho_i = bench_rho(i);
        theta_i = bench_theta(i);
        l = (i == 1) * 2.86 + (i > 1) * 1.65;
        % rho_iplus1 = x(1)
        % theta_iplus1 = x(2)
        rho_iplus1_of_i = @(x) ...
            [ ...
             x(1) - rho_i - k * (x(2) - theta_i); ...
             rho_i ^ 2 + x(1) ^ 2 - 2 * rho_i * x(1) * cos(x(2) - theta_i) - l ^ 2];
        x0 = [rho_i + k * (l / rho_i), theta_i + (l / rho_i)];
        sol = fsolve(rho_iplus1_of_i, x0, options);
        bench_rho(i + 1) = sol(1);
        bench_theta(i + 1) = sol(2);

    end

    result_ki = zeros(1, bench_numb);
    result_x = zeros(1, bench_numb);
    result_y = zeros(1, bench_numb);
    result_v = zeros(1, bench_numb);
    result_v(1) = 1;

    result_alpha = zeros(1, bench_numb - 1);
    result_beta = zeros(1, bench_numb - 1);
    result_k = zeros(1, bench_numb - 1);
    result_b = zeros(1, bench_numb - 1);

    for i = 1:bench_numb
        th = bench_theta(i);
        result_ki(i) = ((r0 + k * th) * cos(th) + k * sin(th)) / (k * cos(th) - (k * th + r0) * sin(th));
        result_x(i) = bench_rho(i) * cos(th);
        result_y(i) = bench_rho(i) * sin(th);
    end

    for i = 1:bench_numb - 1
        result_k(i) = (result_y(i + 1) - result_y(i)) / (result_x(i + 1) - result_x(i));
    end

    for i = 1:bench_numb - 1
        result_alpha(i) = atan(abs((result_ki(i) - result_k(i)) / (result_ki(i) * result_k(i) +1)));
        result_beta(i) = atan(abs((result_ki(i + 1) - result_k(i)) / (result_ki(i + 1) * result_k(i) +1)));
    end

    for i = 1:bench_numb - 1
        result_v(i + 1) = result_v(i) * cos(result_alpha(i)) / cos(result_beta(i));
    end

    for i = 1:bench_numb - 1
        result_b(i) = result_y(i) - result_k(i) * result_x(i);

    end

    result_k1 = (0.275 * result_k(1) + 0.15) / (0.275 - 0.15 * result_k(1));
    result_k2 = (0.275 * result_k(1) - 0.15) / (0.275 + 0.15 * result_k(1));

    % % * 计算x1, y1, x2, y2
    % fun_x1_y1 = @(x1, y1, k, b, k1, xb, yb) ...
    %     [ ...
    %      y1 - k1 * (x1 - xb) - yb; ...
    %      (k * x1 + b - y1) ^ 2 - (0.15 * 0.15 * (k ^ 2 + 1))];

    % options = optimoptions( ...
    %     'fsolve', ...
    %     'Display', 'off', ...
    %     'MaxIterations', 1000, ...
    %     'MaxFunctionEvaluations', 1000, ...
    %     'FunctionTolerance', 1e-10, ...
    %     'StepTolerance', 1e-10);

    % fun_solve_1 = @(xy) fun_x1_y1(xy(1), xy(2), result_k(1), result_b(1), result_k1(1), result_x(1), result_y(1));

    % % init_x1_y1 = [result_x(1), result_y(1); result_x(1) + 0.1, result_y(1) + 0.1];
    % x1_y1 = fsolve(fun_solve_1, init_x1_y1, options);
    % result_x1 = x1_y1(1); % x1
    % result_y1 = x1_y1(2); % y1

    % fun_solve_2 = @(xy) fun_x1_y1(xy(1), xy(2), result_k(2), result_b(1), result_k2(1), result_x(1), result_y(1));

    % % init_x2_y2 = [result_x(2), result_y(2); result_x(2) + 0.1, result_y(2) + 0.1];
    % x2_y2 = fsolve(fun_solve_2, init_x2_y2, options);
    % result_x2 = x2_y2(1); % x2
    % result_y2 = x2_y2(2); % y2

    %* 方程一 y1 = k1x1 + (yb - k1xb)
    %* 方程二 y1 = kx1 + (b +- 0.15 sqrt(k^2 + 1))
    fun1_k = result_k1;
    fun1_b = result_y(1) - result_k1 * result_x(1);
    fun2_k = result_k(1);
    fun2_b = result_b(1) + 0.15 * sqrt(result_k(1) ^ 2 +1);

    result_y1_1 = (fun2_k * fun1_b - fun1_k * fun2_b) / (fun2_k - fun1_k);
    result_x1_1 = (result_y1_1 - fun1_b) / fun1_k;

    fun2_b = result_b(1) - 0.15 * sqrt(result_k(1) ^ 2 +1);

    result_y1_2 = (fun2_k * fun1_b - fun1_k * fun2_b) / (fun2_k - fun1_k);
    result_x1_2 = (result_y1_2 - fun1_b) / fun1_k;

    if (result_x1_1^2 + result_y1_1^2) > (result_x1_2^2 + result_y1_2^2)
        result_x1 = result_x1_1;
        result_y1 = result_y1_1;

    else
        result_x1 = result_x1_2;
        result_y1 = result_y1_2;

    end

    fun1_k = result_k2;
    fun1_b = result_y(1) - result_k2 * result_x(1);
    fun2_k = result_k(1);
    fun2_b = result_b(1) + 0.15 * sqrt(result_k(1) ^ 2 +1);

    result_y2_1 = (fun2_k * fun1_b - fun1_k * fun2_b) / (fun2_k - fun1_k);
    result_x2_1 = (result_y2_1 - fun1_b) / fun1_k;

    fun2_b = result_b(1) - 0.15 * sqrt(result_k(1) ^ 2 +1);

    result_y2_2 = (fun2_k * fun1_b - fun1_k * fun2_b) / (fun2_k - fun1_k);
    result_x2_2 = (result_y2_2 - fun1_b) / fun1_k;

    if (result_x2_1^2 + result_y2_1^2) > (result_x2_2^2 + result_y2_2^2)
        result_x2 = result_x2_1;
        result_y2 = result_y2_1;
    else
        result_x2 = result_x2_2;
        result_y2 = result_y2_2;

    end

    result_d = zeros(bench_numb - 3, 2);

    for i = 1:bench_numb - 3
        result_d(i, 1) = abs((result_k(i + 2) * result_x1 - result_y1 + result_b(i + 2))) / sqrt(result_k(i + 2) ^ 2 +1);
        result_d(i, 2) = abs((result_k(i + 2) * result_x2 - result_y2 + result_b(i + 2))) / sqrt(result_k(i + 2) ^ 2 +1);
    end

    result_d_all(:, :, binary_times) = result_d;
    min_d = min(min(result_d));
    disp(min_d);

    if min_d < 0.15
        disp('NO');
        t_last = t_mid;
    else
        disp('YES');
        t_first = t_mid;
    end

    disp(t_mid);

end
