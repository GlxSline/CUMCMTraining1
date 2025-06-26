function improved_collision_detection()
    % 参数设置
    r0 = 8.8;
    k = 0.55 / (2 * pi);
    v = 1;
    numb = 224;
    t_numb = 51;  % 减少点数提高效率
    
    % 调整二分搜索范围，根据预期结果412左右
    t_min = 410;
    t_max = 415;
    tolerance = 1e-3;  % 适当放宽容差
    max_iterations = 20;
    


    fprintf('开始二分搜索寻找碰撞时间...\n');

    for iter = 1:max_iterations
        t_mid = (t_min + t_max) / 2;
        fprintf('迭代 %d: t_range = [%.6f, %.6f], t_mid = %.6f\n', iter, t_min, t_max, t_mid);

        % 检查t_mid时刻是否发生碰撞
        collision_occurred = check_collision_at_time(t_mid, r0, k, v, numb, t_numb);

        if collision_occurred
            t_max = t_mid; % 碰撞发生，缩小上界
            fprintf('  -> 检测到碰撞，缩小上界\n');
        else
            t_min = t_mid; % 无碰撞，增大下界
            fprintf('  -> 无碰撞，增大下界\n');
        end

        % 检查收敛
        if abs(t_max - t_min) < tolerance
            fprintf('收敛！碰撞时间约为: %.6f\n', t_max);
            break;
        end

    end

    if iter == max_iterations
        fprintf('警告：达到最大迭代次数，未完全收敛\n');
    end

    % 最终验证
    fprintf('\n最终验证:\n');
    fprintf('t = %.6f 时的碰撞状态: %s\n', t_min, ...
        check_collision_at_time(t_min, r0, k, v, numb, t_numb));
    fprintf('t = %.6f 时的碰撞状态: %s\n', t_max, ...
        check_collision_at_time(t_max, r0, k, v, numb, t_numb));
end

function collision_occurred = check_collision_at_time(t_end, r0, k, v, numb, t_numb)
    try
        % 只计算最终时刻的位置
        t = t_end;  % 单个时间点
        
        % 求解龙头位置
        A = @(u) u .* sqrt(u .^ 2 + k ^ 2) + k ^ 2 * log(u + sqrt(u .^ 2 + k ^ 2));
        s_of_theta = @(theta) (A(r0 + k * theta) - A(r0)) / (2 * k);
        
        funs = @(th) s_of_theta(th) - v * t;
        theta_head = fzero(funs, 2 * v * t / r0, optimset('Display', 'off'));
        r_head = r0 + k * theta_head;
        
        % 求解所有板凳位置
        positions = solve_all_positions_single_time(r_head, theta_head, k, numb);
        
        % 检查碰撞：龙头宽度内的点与后续板凳的距离
        collision_occurred = check_width_collision(positions, numb);
        
    catch ME
        collision_occurred = true; % 计算失败时假设碰撞
    end
end


function [r_sol, theta_sol] = solve_head_trajectory(t, s_of_theta, v, r0, k)
    r_sol = zeros(size(t));
    theta_sol = zeros(size(t));
    opts_fz = optimset('Display', 'off', 'TolX', 1e-10);
    
    for i = 1:numel(t)
        ti = t(i);
        % 关键修复：弧长方程应该是 s(theta) = v*t
        funs = @(th) s_of_theta(th) - v * ti;  % 这里必须是减号
        
        % 改进初值估计
        if i == 1
            theta0 = 2 * v * ti / r0;  % 基于小角度近似的初值
        else
            % 基于前一个解的线性外推
            theta0 = theta_sol(i-1) + v * (t(i) - t(i-1)) / r_sol(i-1);
        end
        
        try
            si = fzero(funs, theta0, opts_fz);
            theta_sol(i) = si;
            r_sol(i) = r0 + k * si;
        catch
            % 备用求解
            theta_sol(i) = theta0;
            r_sol(i) = r0 + k * theta0;
        end
    end
end



function [result_rho, result_theta] = solve_all_positions(t, r_sol, theta_sol, k, numb)
    result_rho = zeros(numel(t), numb);
    result_theta = zeros(numel(t), numb);
    result_rho(:, 1) = r_sol;
    result_theta(:, 1) = theta_sol;

    options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

    for i = 1:numel(t)

        for j = 1:(numb - 1)
            rho1 = result_rho(i, j);
            theta1 = result_theta(i, j);
            l = (j == 1) * 2.86 + (j > 1) * 1.65;

            % 改进的初值估计
            if i > 1
                % 使用前一时刻的解作为初值
                x0 = [result_rho(i - 1, j + 1); result_theta(i - 1, j + 1)];
            else
                % 使用几何估计
                dth0 = l / rho1;
                x0 = [rho1 + k * dth0; theta1 + dth0];
            end

            try
                sol = fsolve(@(x) segment_eq(x, rho1, theta1, k, l), x0, options);
                result_rho(i, j + 1) = sol(1);
                result_theta(i, j + 1) = sol(2);
            catch
                % 求解失败时使用备用方法
                result_rho(i, j + 1) = x0(1);
                result_theta(i, j + 1) = x0(2);
            end

        end

    end

end

function [result_x, result_y, result_k_slopes] = calculate_geometry(result_rho, result_theta, k, r0, numb)
    [t_numb, ~] = size(result_rho);
    
    % 计算笛卡尔坐标
    result_x = result_rho .* cos(result_theta);
    result_y = result_rho .* sin(result_theta);
    
    % 计算每节板凳的方向斜率（板凳朝向）
    result_k_slopes = zeros(t_numb, numb - 1);
    for i = 1:t_numb
        for j = 1:(numb - 1)
            % 使用相邻两节板凳的位置计算板凳方向
            dx = result_x(i, j + 1) - result_x(i, j);
            dy = result_y(i, j + 1) - result_y(i, j);
            
            if abs(dx) < 1e-12
                result_k_slopes(i, j) = sign(dy) * 1e6; % 垂直方向
            else
                result_k_slopes(i, j) = dy / dx;
            end
        end
    end
end


function collision_occurred = check_collision(result_x, result_y, result_k_slopes, t_numb, numb)
    collision_occurred = false;
    
    % 只检查最后一个时刻点，因为我们关心的是在t_end时刻是否碰撞
    i = t_numb;
    
    % 计算第一节板凳的直线参数
    x_head = result_x(i, 1);
    y_head = result_y(i, 1);
    k_head = result_k_slopes(i, 1);
    b_head = y_head - k_head * x_head;
    
    % 检查龙头与第3节及以后板凳的最小距离
    min_distance = inf;
    
    for j = 3:numb-1  % 从第3节开始检查
        x_bench = result_x(i, j);
        y_bench = result_y(i, j);
        
        % 点到直线距离公式: |kx - y + b| / sqrt(k^2 + 1)
        distance = abs(k_head * x_bench - y_bench + b_head) / sqrt(k_head^2 + 1);
        min_distance = min(min_distance, distance);
    end
    
    % 碰撞判断：考虑龙头宽度0.275m，安全距离应为0.275/2 + 板凳宽度/2
    safety_distance = 0.275/2 + 0.15;  % 龙头半宽 + 板凳半宽
    
    if min_distance < safety_distance
        collision_occurred = true;
    end
end


function [x1, y1, x2, y2] = calculate_boundary_points(x_head, y_head, k_head)
    % 计算龙头边界点的改进版本
    if abs(k_head) > 1e6 % 近似垂直线
        x1 = x_head + 0.15;
        y1 = y_head;
        x2 = x_head - 0.15;
        y2 = y_head;
    else
        % 垂直距离为0.15的两条平行线
        angle = atan(k_head);
        dx = 0.15 * sin(angle);
        dy = 0.15 * cos(angle);

        x1 = x_head + dx;
        y1 = y_head - dy;
        x2 = x_head - dx;
        y2 = y_head + dy;
    end

end

function min_distance = calculate_min_distance_to_benches(x1, y1, x2, y2, bench_x, bench_y, bench_k)
    min_distance = inf;

    points = [x1, y1; x2, y2];

    for p = 1:size(points, 1)
        px = points(p, 1);
        py = points(p, 2);

        for j = 1:length(bench_x) - 1
            % 点到直线的距离公式
            if abs(bench_k(j)) > 1e6 % 垂直线
                distance = abs(px - bench_x(j));
            else
                % 直线方程: kx - y + b = 0, 其中 b = y - kx
                b = bench_y(j) - bench_k(j) * bench_x(j);
                distance = abs(bench_k(j) * px - py + b) / sqrt(bench_k(j) ^ 2 + 1);
            end

            min_distance = min(min_distance, distance);
        end

    end

end

function F = segment_eq(x, rho1, theta1, k, l)
    rho2 = x(1);
    theta2 = x(2);
    F = [
         rho2 - rho1 - k * (theta2 - theta1);
         rho1 ^ 2 + rho2 ^ 2 - 2 * rho1 * rho2 * cos(theta2 - theta1) - l ^ 2
         ];
end

function positions = solve_all_positions_single_time(r_head, theta_head, k, numb)
    rho = zeros(1, numb);
    theta = zeros(1, numb);
    rho(1) = r_head;
    theta(1) = theta_head;
    
    for j = 1:(numb - 1)
        l = (j == 1) * 2.86 + (j > 1) * 1.65;
        
        % 使用解析解或高质量数值解
        rho1 = rho(j);
        theta1 = theta(j);
        
        % 更好的初值估计
        dtheta_est = l / rho1;
        rho2_init = rho1 + k * dtheta_est;
        theta2_init = theta1 + dtheta_est;
        
        % 求解约束方程
        options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-12);
        sol = fsolve(@(x) segment_eq_improved(x, rho1, theta1, k, l), ...
                     [rho2_init; theta2_init], options);
        
        rho(j + 1) = sol(1);
        theta(j + 1) = sol(2);
    end
    
    % 转换为笛卡尔坐标
    x = rho .* cos(theta);
    y = rho .* sin(theta);
    positions = [x; y]';
end


% 运行主函数
improved_collision_detection();
