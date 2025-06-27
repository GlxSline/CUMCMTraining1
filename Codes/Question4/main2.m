%* 龙头坐标方程

% function F = theta2_of_theta1(input)

% end

%* 约束方程

function F = segment_eq_1(x, rho1, theta1, k, l)
    rho2 = x(1);
    theta2 = x(2);
    F = [
         rho2 - rho1 - k * (theta2 - theta1);
         rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
         ];

end

% function F = segment_eq_2(x, rho1, theta1, R, xE, yE, l)
%     rho2 = x(1);
%     theta2 = x(2);
%     x1 = rho1 * cos (theta1);
%     x2 = rho2 * cos (theta2);
%     y1 = rho1 * sin (theta1);
%     y2 = rho2 * sin (theta2);
%     % 避免除零错误
%     if abs(x1 - xE) < 1e-10 || abs(x2 - xE) < 1e-10
%         F = [0; 0];
%         return;
%     end

%     k1 = (y1 - yE) / (x1 - xE);
%     k2 = (y2 - yE) / (x2 - xE);

%     vector1 = [x1 - xE; y1 - yE];
%     vector2 = [x2 - xE; y2 - yE];
%     % l1 = sqrt((y1 - yE)^2 + (x1 - xE)^2);
%     % l2 = sqrt((y2 - yE)^2 + (x2 - xE)^2);
%     if abs(1 + k1 * k2) < 1e-10
%         theta = pi / 2; % 垂直情况
%     else
%         theta = atan((k1 - k2) / (1 + k1 * k2));
%     end

%     % theta = atan((k1 - k2) / (1 + k1 * k2));

%     rotation_matrix = [cos(theta), sin(theta); -1 * sin(theta), cos(theta)];

%     F = vector2 - rotation_matrix * vector1;
%     %  R ^ 2 + R ^ 2 - 2 * R * R * cos(theta) - l ^ 2;

% end

% function F = segment_eq_3(x, rho1, theta1, R, xE, yE, l)
%     rho2 = x(1);
%     theta2 = x(2);
%     x1 = rho1 * cos (theta1);
%     x2 = rho2 * cos (theta2);
%     y1 = rho1 * sin (theta1);
%     y2 = rho2 * sin (theta2);

%     if abs(x1 - xE) < 1e-10 || abs(x2 - xE) < 1e-10
%         F = [0; 0];
%         return;
%     end

%     k1 = (y1 - yE) / (x1 - xE);
%     k2 = (y2 - yE) / (x2 - xE);

%     vector1 = [x1 - xE; y1 - yE];
%     vector2 = [x2 - xE; y2 - yE];
%     % l1 = sqrt((y1 - yE)^2 + (x1 - xE)^2);
%     % l2 = sqrt((y2 - yE)^2 + (x2 - xE)^2);
%     if abs(1 + k1 * k2) < 1e-10
%         theta = pi / 2; % 垂直情况
%     else
%         theta = atan((k1 - k2) / (1 + k1 * k2));
%     end

%     % theta = atan((k1 - k2) / (1 + k1 * k2));

%     rotation_matrix = [cos(theta), -1 * sin(theta); sin(theta), cos(theta)];

%     F = vector2 - rotation_matrix * vector1;
%     %  R ^ 2 + R ^ 2 - 2 * R * R * cos(theta) - l ^ 2;

% end

function F = segment_eq_4(x, rho1, theta1, k, l)
    %* ρ = k (θ + pi)
    rho2 = x(1);
    theta2 = x(2);
    F = [
        %  rho2 - rho1 - k * (theta2 - theta1);
        rho2 - k * (theta2 + pi);
         rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
         ];

end

function F = segment_eq_1_2(x, rho1, theta1, k, l)

    rho2 = x(1);
    theta2 = x(2);
    x2 = rho2 * cos(theta2);
    y2 = rho2 * sin(theta2);
    x1 = rho1 * cos(theta1);
    y1 = rho1 * sin(theta1);

    F = [
         (x2 - x1) ^ 2 + (y2 - y1) ^ 2 - l ^ 2;
         rho2 - k * theta2;
         ];

end

function F = segment_eq_2_3(x, x1, y1, x_E2, y_E2, R, l)
    x2 = x(1);
    y2 = x(2);

    F = [
         (x2 - x1) ^ 2 + (y2 - y1) ^ 2 - l ^ 2;
         (x2 - x_E2) ^ 2 + (y2 - y_E2) ^ 2 - R ^ 2;
         ];

end

function F = segment_eq_3_4(x, x1, y1, x_E4, y_E4, R, l)
    x2 = x(1);
    y2 = x(2);

    F = [
         (x2 - x1) ^ 2 + (y2 - y1) ^ 2 - l ^ 2;
         (x2 - x_E4) ^ 2 + (y2 - y_E4) ^ 2 - R ^ 2;
         ];

end

%* 计算x_E, y_E
theta_E = 4.5 * 2 * pi / 1.7;
rho_E = 4.5;
x_E1 = rho_E * cos(theta_E);
y_E1 = rho_E * sin(theta_E);
x_E3 = x_E1 / (-3);
y_E3 = y_E1 / (-3);
x_E5 = -1 * x_E1;
y_E5 = -1 * y_E1;

k_E1 = (sin(theta_E) + theta_E * cos(theta_E)) / (cos(theta_E) - theta_E * sin(theta_E));
k_E1E5 = y_E1 / x_E1;
k_E1E2 = -1 / k_E1;
theta_E2 = atan(k_E1E2);

gamma_E = atan((k_E1E5 - k_E1E2) / (1 + k_E1E2 * k_E1E5));
r_E1E2 = 3 / cos(gamma_E);
r_E3E4 = r_E1E2 / 2;
r_E1E3 = sqrt((x_E1 - x_E3) ^ 2 + (y_E1 - y_E3) ^ 2);

x_E2 = x_E1 + r_E1E2 * cos(theta_E2);
y_E2 = y_E1 + r_E1E2 * sin(theta_E2);
x_E4 = (3 * x_E3 - x_E2) / 2;
y_E4 = (3 * y_E3 - y_E2) / 2;

k_E2E3 = (y_E3 - y_E2) / (x_E3 - x_E2);

width = 1.7;
k = width / (2 * pi);
v = 1;
bench_numb = 224;
t = linspace(0, 100, 101);
r0 = 4.5;
theta0 = r0 / k - pi;

%* t1
theta_E1E3 = pi + atan((k_E1E2 - k_E2E3) / (1 + k_E1E2 * k_E2E3));
l_E1E3 = theta_E1E3 * r_E1E2;
t1 = l_E1E3 / v;

%* t2
l_E3E5 = theta_E1E3 * r_E3E4;
t2 = t1 + l_E3E5 / v;

%* x1, y1, rho1, theta1
result_rho = zeros(numel(t), bench_numb);
result_theta = zeros(numel(t), bench_numb);
result_x = zeros(numel(t), bench_numb);
result_y = zeros(numel(t), bench_numb);

for i = 1:numel(t)
    % x = rho * cos (theta);
    % y = rho * sin (theta);
    % r_AE1 = sqrt((x - x_E1) ^ 2 + (y - y_E1) ^ 2);
    % k_AE1 = (y - y_E1) / (x - x_E1);

    ta = t(i);

    if ta <= t1 %* II
        theta_AE1 = v * ta / r_E1E2;
        vector_E2E1 = [x_E1 - x_E2; y_E1 - y_E2];
        rotation_matrix = [cos(theta_AE1), sin(theta_AE1); -1 * sin(theta_AE1), cos(theta_AE1)];
        vector_E2A = rotation_matrix * vector_E2E1;
        x_i = vector_E2A(1) + x_E2;
        y_i = vector_E2A(2) + y_E2;
        rho_i = sqrt(x_i ^ 2 + y_i ^ 2);
        theta_i = atan2(y_i, x_i);

        % elseif (ta > t1 && t(i - 1) < t1)

    elseif ta <= t2 %* III
        theta_AE3 = v * (ta - t1) / r_E3E4;
        vector_E4E3 = [x_E3 - x_E4; y_E3 - y_E4];
        rotation_matrix = [cos(theta_AE3), -1 * sin(theta_AE3); sin(theta_AE3), cos(theta_AE3)];
        vector_E4A = rotation_matrix * vector_E4E3;
        x_i = vector_E4A(1) + x_E4;
        y_i = vector_E4A(2) + y_E4;
        rho_i = sqrt(x_i ^ 2 + y_i ^ 2);
        theta_i = atan2(y_i, x_i);

        % elseif (ta > t2 && t(i - 1) < t2)

    else %* IV
        opts_fz = optimset('Display', 'off');
        theta_start = result_theta(i - 1, 1);
        % t_of_theta2 = @(theta_i) (k / 2) * ((theta_i+ pi) * sqrt(1 + (theta_i+ pi) ^ 2) - pi * sqrt(1 + pi ^ 2) ...
        %     + log((theta_i+ pi) + sqrt(1 + (theta_i+ pi) ^ 2)) - log(pi + sqrt(1 + pi ^ 2))) - v * (ta - t2);

        t_of_theta2 = @(theta_i) (k / 2) * ...
            ((theta_i + pi) * sqrt(1 + (theta_i + pi) ^ 2) - ...
            (theta0 + pi) * sqrt(1 + (theta0 + pi) ^ 2) + ...
            log((theta_i + pi) + sqrt(1 + (theta_i + pi) ^ 2)) - ...
            log((theta0 + pi) + sqrt(1 + (theta0 + pi) ^ 2))) - v * (ta - t2);

        theta_i = fzero(t_of_theta2, theta_start, opts_fz);
        rho_i = k * (theta_i + pi);
        x_i = rho_i * cos(theta_i);
        y_i = rho_i * sin(theta_i);

    end

    % for j = i:numel(t)
    %     %* r = k (θ + Π/2)；
    %     opts_fz = optimset('Display', 'off');
    %     t_of_theta2 = @(theta_i) (k / 2) * (theta_i * sqrt(1 + theta_i ^ 2) - pi * sqrt(1 + pi ^ 2) ...
    %         + log(theta_i + sqrt(1 + theta_i ^ 2)) - log(pi + sqrt(1 + pi ^ 2))) - v * ta;

    %     initial_guess = pi / 2;
    %     theta_i = fzero(t_of_theta2, pi / 2, opts_fz);
    %     rho_i = k * (theta_i + pi / 2);
    %     x_i = rho_i * cos(theta_i);
    %     y_i = rho_i * sin(theta_i);
    % end

    result_rho(i, 1) = rho_i;
    result_theta(i, 1) = theta_i;
    result_x(i, 1) = x_i;
    result_y(i, 1) = y_i;

end

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

for i = 1:numel(t)

    for j = 1:bench_numb - 1
        rho1 = result_rho(i, j);
        theta1 = result_theta(i, j);
        l = (j == 1) * 2.86 + (j > 1) * 1.65;

        area_numb = calculate_area(rho1, theta1, x_E2, y_E2, r_E1E2, k_E1E2, x_E4, y_E4, r_E3E4, k_E2E3, x_E5, y_E5, k, l);
        options = optimoptions('fsolve', 'Display', 'off');

        if (area_numb(1) == 1 && area_numb(2) == 1)
            rho1 = result_rho(i, j);
            theta1 = result_theta(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;
            x0 = [result_rho(i, j), result_theta(i, j)];
            sol = fsolve(@(x) segment_eq_1(x, rho1, theta1, k, l), x0, options);
            result_rho(i, j + 1) = sol(1);
            result_theta(i, j + 1) = sol(2);
            result_x(i, j + 1) = sol(1) * cos(sol(2));
            result_y(i, j + 1) = sol(1) * sin(sol(2));

        elseif (area_numb(1) == 2 && area_numb(2) == 2)
            rho1 = result_rho(i, j);
            theta1 = result_theta(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;
            % rotation_matrix_temp = [cos(0.9), sin(0.9); -1 * sin(0.9), cos(0.9)];
            % x0 = rotation_matrix_temp * [result_rho(i, j); result_theta(i, j)];
            % sol = fsolve(@(x) segment_eq_2(x, rho1, theta1, r_E1E2, x_E2, y_E2, l), x0, options);
            % result_rho(i, j + 1) = sol(1);
            % result_theta(i, j + 1) = sol(2);
            % result_x(i, j + 1) = sol(1) * cos(sol(2));
            % result_y(i, j + 1) = sol(1) * sin(sol(2));

            theta_rotation = acos((2 * r_E1E2 * r_E1E2 - l ^ 2) / (2 * r_E1E2 * r_E1E2));

            x_1 = rho1 * cos (theta1);
            y_1 = rho1 * sin (theta1);

            vector1 = [x_1 - x_E2; y_1 - y_E2];
            rotation_matrix = [cos(theta_rotation), -1 * sin(theta_rotation); sin(theta_rotation), cos(theta_rotation)];
            vector2 = rotation_matrix * vector1;
            result_x(i, j + 1) = vector2(1) + x_E2;
            result_y(i, j + 1) = vector2(2) + y_E2;
            result_rho(i, j + 1) = sqrt(result_x(i, j + 1) ^ 2 + result_y(i, j + 1) ^ 2);
            result_theta(i, j + 1) = atan2(result_y(i, j + 1), result_x(i, j + 1));

        elseif (area_numb(1) == 3 && area_numb(2) == 3)
            rho1 = result_rho(i, j);
            theta1 = result_theta(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;
            % x0 = [result_rho(i, j), result_theta(i, j)];
            % sol = fsolve(@(x) segment_eq_3(x, rho1, theta1, r_E3E4, x_E4, y_E4, l), x0, options);
            % result_rho(i, j + 1) = sol(1);
            % result_theta(i, j + 1) = sol(2);
            % result_x(i, j + 1) = sol(1) * cos(sol(2));
            % result_y(i, j + 1) = sol(1) * sin(sol(2));
            theta_rotation = acos((2 * r_E3E4 * r_E3E4 - l ^ 2) / (2 * r_E3E4 * r_E3E4));

            x_1 = rho1 * cos (theta1);
            y_1 = rho1 * sin (theta1);

            vector1 = [x_1 - x_E4; y_1 - y_E4];
            rotation_matrix = [cos(theta_rotation), sin(theta_rotation); -1 * sin(theta_rotation), cos(theta_rotation)];
            vector2 = rotation_matrix * vector1;
            result_x(i, j + 1) = vector2(1) + x_E4;
            result_y(i, j + 1) = vector2(2) + y_E4;
            result_rho(i, j + 1) = sqrt(result_x(i, j + 1) ^ 2 + result_y(i, j + 1) ^ 2);
            result_theta(i, j + 1) = atan2(result_y(i, j + 1), result_x(i, j + 1));

        elseif (area_numb(1) == 4 && area_numb(2) == 4)
            rho1 = result_rho(i, j);
            theta1 = result_theta(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;
            theta0 = l / result_rho(i, j);
            % rotation_matrix = [cos(theta0), sin(theta0); -1 * sin(theta0), cos(theta0)];
            % x_rotation = rotation_matrix * [result_x(i, j); result_y(i, j)];
            % x0(1) = sqrt(x_rotation(1) ^ 2 + x_rotation(2) ^ 2);
            % x0(2) = x0(1) / k - pi;
            x0 = [rho1, theta1 - theta0];
            sol = fsolve(@(x) segment_eq_4(x, rho1, theta1, k, l), x0, options);
            result_rho(i, j + 1) = sol(1);
            result_theta(i, j + 1) = sol(2);
            result_x(i, j + 1) = sol(1) * cos(sol(2));
            result_y(i, j + 1) = sol(1) * sin(sol(2));

            % fprintf('4-4');
            % disp(i);disp(j);

        elseif (area_numb(1) == 2 && area_numb(2) == 1)
            rho1 = result_rho(i, j);
            theta1 = result_theta(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;

            if i == 1
                x0 = [4.67, 17.30];
            else
                x0 = [result_rho(i - 1, j + 1), result_theta(i - 1, j + 1)];
            end

            % fprintf('1-2');
            % disp(i);disp(j);

            sol = fsolve(@(x) segment_eq_1_2(x, rho1, theta1, k, l), x0, options);
            result_rho(i, j + 1) = sol(1);
            result_theta(i, j + 1) = sol(2);
            result_x(i, j + 1) = sol(1) * cos(sol(2));
            result_y(i, j + 1) = sol(1) * sin(sol(2));

        elseif (area_numb(1) == 3 && area_numb(2) == 2)
            x1 = result_x(i, j);
            y1 = result_y(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;
            % x0 = [x1, y1];
            x0 = [0.488, 1.428];
            sol = fsolve(@(x) segment_eq_2_3(x, x1, y1, x_E2, y_E2, r_E1E2, l), x0, options);
            result_x(i, j + 1) = sol(1);
            result_y(i, j + 1) = sol(2);
            result_rho(i, j + 1) = sqrt(result_x(i, j + 1) ^ 2 + result_y(i, j + 1) ^ 2);
            result_theta(i, j + 1) = atan2(result_y(i, j + 1), result_x(i, j + 1));

        elseif (area_numb(1) == 4 && area_numb(2) == 3)
            x1 = result_x(i, j);
            y1 = result_y(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;
            % x0 = [result_x(i - 1, j + 1), result_y(i - 1, j + 1)];
            x0 = [3.049, 1.718];
            sol = fsolve(@(x) segment_eq_3_4(x, x1, y1, x_E4, y_E4, r_E3E4, l), x0, options);
            result_x(i, j + 1) = sol(1);
            result_y(i, j + 1) = sol(2);
            result_rho(i, j + 1) = sqrt(result_x(i, j + 1) ^ 2 + result_y(i, j + 1) ^ 2);
            result_theta(i, j + 1) = atan2(result_y(i, j + 1), result_x(i, j + 1));

            % fprintf('3-4');
            % disp(i);disp(j);

        end

    end

    % 生成动画
    % animate_dragon_motion(t, result_x, result_y, bench_numb, ...
    %     x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, ...
    %     r_E1E2, r_E3E4);

end

% 在你的主循环结束后添加：
% 绘制特定时间步

time_indices = linspace(1, 30, 30); % 选择要显示的时间步

for idx = time_indices

    if idx <= length(t)
        plot_dragon_at_time(idx, t, result_x, result_y, bench_numb, ...
            x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, ...
            r_E1E2, r_E3E4);
        pause(0.5); % 停留2秒观察
    end

end
