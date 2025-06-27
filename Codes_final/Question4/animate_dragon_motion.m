function animate_dragon_motion(t, result_x, result_y, bench_numb, x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, r_E1E2, r_E3E4)
    % 创建动画显示龙的运动过程

    figure('Position', [100, 100, 1200, 800]);

    % 预先计算坐标轴范围
    all_x = result_x(result_x ~= 0);
    all_y = result_y(result_y ~= 0);
    special_x = [x_E1, x_E2, x_E3, x_E4, x_E5];
    special_y = [y_E1, y_E2, y_E3, y_E4, y_E5];
    margin = 3;
    x_range = [min([all_x(:); special_x(:)]) - margin, max([all_x(:); special_x(:)]) + margin];
    y_range = [min([all_y(:); special_y(:)]) - margin, max([all_y(:); special_y(:)]) + margin];

    for i = 1:5:length(t) % 每5步显示一次，加快动画速度
        clf;
        hold on;
        grid on;
        axis equal;

        % 绘制轨道和特殊点（静态背景）
        plot_background(x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, r_E1E2, r_E3E4);

        % 绘制龙身
        plot_dragon_body(i, result_x, result_y, bench_numb);

        % 设置图形属性
        xlim(x_range);
        ylim(y_range);
        xlabel('X坐标', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Y坐标', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('舞龙动画 - 时间: t = %.2f s', t(i)), ...
            'FontSize', 14, 'FontWeight', 'bold');

        hold off;
        pause(0.1); % 控制动画速度
        drawnow;
    end

end

function plot_background(x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, r_E1E2, r_E3E4)
    % 绘制静态背景元素

    % 螺旋轨道
    theta_spiral = linspace(-pi, 6 * pi, 1000);
    k = 1.7 / (2 * pi);
    rho_spiral = k * (theta_spiral + pi);
    x_spiral = rho_spiral .* cos(theta_spiral);
    y_spiral = rho_spiral .* sin(theta_spiral);
    plot(x_spiral, y_spiral, 'k--', 'LineWidth', 1, 'Color', [0.8 0.8 0.8]);

    % 圆弧轨道
    theta_arc = linspace(0, 2 * pi, 100);
    x_arc1 = x_E2 + r_E1E2 * cos(theta_arc);
    y_arc1 = y_E2 + r_E1E2 * sin(theta_arc);
    plot(x_arc1, y_arc1, '--', 'LineWidth', 1, 'Color', [0.7 0.7 1]);

    x_arc2 = x_E4 + r_E3E4 * cos(theta_arc);
    y_arc2 = y_E4 + r_E3E4 * sin(theta_arc);
    plot(x_arc2, y_arc2, '--', 'LineWidth', 1, 'Color', [0.7 1 0.7]);

    % 特殊点
    special_x = [x_E1, x_E2, x_E3, x_E4, x_E5];
    special_y = [y_E1, y_E2, y_E3, y_E4, y_E5];
    plot(special_x, special_y, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'red');
end

function plot_dragon_body(i, result_x, result_y, bench_numb)
    % 绘制龙身

    % 获取有效节点
    valid_indices = (result_x(i, :) ~= 0) & (result_y(i, :) ~= 0);
    valid_nodes = find(valid_indices);

    if length(valid_nodes) < 2
        return;
    end

    % 绘制连线（带渐变颜色效果）
    n_valid = length(valid_nodes);

    for j = 1:(n_valid - 1)
        idx1 = valid_nodes(j);
        idx2 = valid_nodes(j + 1);

        % 根据位置计算颜色（龙头为红色，龙尾为蓝色）
        color_ratio = (j - 1) / (n_valid - 1);
        line_color = [1 - color_ratio, 0, color_ratio];

        plot([result_x(i, idx1), result_x(i, idx2)], ...
            [result_y(i, idx1), result_y(i, idx2)], ...
            '-', 'LineWidth', 2.5, 'Color', line_color);
    end

    % 绘制节点
    plot(result_x(i, valid_indices), result_y(i, valid_indices), ...
        'o', 'MarkerSize', 3, 'MarkerFaceColor', [0.2 0.2 0.8], ...
        'MarkerEdgeColor', 'black', 'LineWidth', 0.5);

    % 突出显示龙头
    if valid_indices(1)
        plot(result_x(i, 1), result_y(i, 1), 'o', 'MarkerSize', 12, ...
            'MarkerFaceColor', 'yellow', 'MarkerEdgeColor', 'red', 'LineWidth', 3);
    end

    % 突出显示龙尾
    last_valid = valid_nodes(end);
    plot(result_x(i, last_valid), result_y(i, last_valid), 's', ...
        'MarkerSize', 8, 'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'blue', 'LineWidth', 2);
end
