function plot_dragon_at_time(i, t, result_x, result_y, bench_numb, x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, r_E1E2, r_E3E4)
    % 绘制第i个时间步的龙头位置和所有节点连线
    % 输入参数：
    % i: 时间步索引
    % t: 时间向量
    % result_x, result_y: 节点坐标矩阵
    % bench_numb: 节点总数
    % x_E*, y_E*: 特殊点坐标
    % r_E1E2, r_E3E4: 圆弧半径
    
    figure('Position', [100, 100, 1200, 800]);
    hold on;
    grid on;
    axis equal;
    
    % 1. 绘制螺旋轨道
    theta_spiral = linspace(-pi, 6*pi, 1000);
    k = 1.7 / (2 * pi);
    rho_spiral = k * (theta_spiral + pi);
    x_spiral = rho_spiral .* cos(theta_spiral);
    y_spiral = rho_spiral .* sin(theta_spiral);
    plot(x_spiral, y_spiral, 'k--', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
    
    % 2. 绘制圆弧轨道
    % 圆弧E1E2 (以E2为圆心)
    theta_arc1 = linspace(0, 2*pi, 100);
    x_arc1 = x_E2 + r_E1E2 * cos(theta_arc1);
    y_arc1 = y_E2 + r_E1E2 * sin(theta_arc1);
    plot(x_arc1, y_arc1, 'b--', 'LineWidth', 1.5, 'Color', [0.5 0.5 1]);
    
    % 圆弧E3E4 (以E4为圆心)
    theta_arc2 = linspace(0, 2*pi, 100);
    x_arc2 = x_E4 + r_E3E4 * cos(theta_arc2);
    y_arc2 = y_E4 + r_E3E4 * sin(theta_arc2);
    plot(x_arc2, y_arc2, 'g--', 'LineWidth', 1.5, 'Color', [0.5 1 0.5]);
    
    % 3. 绘制特殊点E1-E5
    special_points_x = [x_E1, x_E2, x_E3, x_E4, x_E5];
    special_points_y = [y_E1, y_E2, y_E3, y_E4, y_E5];
    plot(special_points_x, special_points_y, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
    
    % 标注特殊点
    text(x_E1, y_E1, '  E1', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E2, y_E2, '  E2', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E3, y_E3, '  E3', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E4, y_E4, '  E4', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E5, y_E5, '  E5', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    
    % 4. 绘制龙身节点连线
    for j = 1:(bench_numb-1)
        if result_x(i, j) ~= 0 && result_y(i, j) ~= 0 && ...
           result_x(i, j+1) ~= 0 && result_y(i, j+1) ~= 0
            % 绘制连线
            plot([result_x(i, j), result_x(i, j+1)], ...
                 [result_y(i, j), result_y(i, j+1)], ...
                 'b-', 'LineWidth', 2);
        end
    end
    
    % 5. 绘制所有节点
    valid_indices = (result_x(i, :) ~= 0) & (result_y(i, :) ~= 0);
    plot(result_x(i, valid_indices), result_y(i, valid_indices), ...
         'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'blue');
    
    % 6. 突出显示龙头（第一个节点）
    if result_x(i, 1) ~= 0 && result_y(i, 1) ~= 0
        plot(result_x(i, 1), result_y(i, 1), 'go', 'MarkerSize', 10, ...
             'MarkerFaceColor', 'green', 'LineWidth', 2);
        text(result_x(i, 1), result_y(i, 1), '  龙头', 'FontSize', 12, ...
             'Color', 'green', 'FontWeight', 'bold');
    end
    
    % 7. 突出显示龙尾（最后一个有效节点）
    last_valid = find(valid_indices, 1, 'last');
    if ~isempty(last_valid)
        plot(result_x(i, last_valid), result_y(i, last_valid), 'mo', ...
             'MarkerSize', 8, 'MarkerFaceColor', 'magenta', 'LineWidth', 2);
        text(result_x(i, last_valid), result_y(i, last_valid), '  龙尾', ...
             'FontSize', 10, 'Color', 'magenta', 'FontWeight', 'bold');
    end
    
    % 8. 设置图形属性
    xlabel('X坐标', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Y坐标', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('舞龙队形图 - 时间: t = %.2f s (步骤 %d/%d)', ...
          t(i), i, length(t)), 'FontSize', 14, 'FontWeight', 'bold');
    
    % 9. 添加图例
    legend({'螺旋轨道', '圆弧轨道1', '圆弧轨道2', '关键点', '龙身连线', ...
            '龙身节点', '龙头', '龙尾'}, ...
           'Location', 'best', 'FontSize', 10);
    
%     % 10. 设置坐标轴范围
%     all_x = [result_x(i, valid_indices), special_points_x];
%     all_y = [result_y(i, valid_indices), special_points_y];
%     if ~isempty(all_x)
%         margin = 2;
%         xlim([min(all_x) - margin, max(all_x) + margin]);
%         ylim([min(all_y) - margin, max(all_y) + margin]);
%     end
    
    hold off;
end
