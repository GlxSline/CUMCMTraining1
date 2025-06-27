function plot_dragon_at_time(i, t, result_x, result_y, bench_numb, x_E1, y_E1, x_E2, y_E2, x_E3, y_E3, x_E4, y_E4, x_E5, y_E5, r_E1E2, r_E3E4)
    
    figure('Position', [100, 100, 1200, 800]);
    hold on;
    grid on;
    axis equal;
    
    theta_spiral = linspace(13.4903684536503, 11*pi, 1000);
    k = 1.7 / (2 * pi);
    rho_spiral = k * (theta_spiral + pi);
    x_spiral = rho_spiral .* cos(theta_spiral);
    y_spiral = rho_spiral .* sin(theta_spiral);
    plot(x_spiral, y_spiral, 'k--', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
 
    theta_spiral = linspace(16.6319611072401, 12*pi, 1000);
    k = 1.7 / (2 * pi);
    rho_spiral = k * (theta_spiral);
    x_spiral = rho_spiral .* cos(theta_spiral);
    y_spiral = rho_spiral .* sin(theta_spiral);
    plot(x_spiral, y_spiral, 'k--', 'LineWidth', 1, 'Color', [0.7 0.7 0.7]);
    
    theta_arc1 = linspace(0, 2*pi, 100);
    x_arc1 = x_E2 + r_E1E2 * cos(theta_arc1);
    y_arc1 = y_E2 + r_E1E2 * sin(theta_arc1);
    plot(x_arc1, y_arc1, 'b--', 'LineWidth', 1, 'Color', [0.5 0.5 1]);
    
    theta_arc2 = linspace(0, 2*pi, 100);
    x_arc2 = x_E4 + r_E3E4 * cos(theta_arc2);
    y_arc2 = y_E4 + r_E3E4 * sin(theta_arc2);
    plot(x_arc2, y_arc2, 'g--', 'LineWidth', 1, 'Color', [0.5 1 0.5]);
    
    special_points_x = [x_E1, x_E2, x_E3, x_E4, x_E5];
    special_points_y = [y_E1, y_E2, y_E3, y_E4, y_E5];
    plot(special_points_x, special_points_y, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'red');
    
    text(x_E1, y_E1, '  E1', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E2, y_E2, '  E2', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E3, y_E3, '  E3', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E4, y_E4, '  E4', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    text(x_E5, y_E5, '  E5', 'FontSize', 10, 'Color', 'red', 'FontWeight', 'bold');
    
    for j = 1:(bench_numb-1)
        if result_x(i, j) ~= 0 && result_y(i, j) ~= 0 && ...
           result_x(i, j+1) ~= 0 && result_y(i, j+1) ~= 0

           plot([result_x(i, j), result_x(i, j+1)], ...
                 [result_y(i, j), result_y(i, j+1)], ...
                 'b-', 'LineWidth', 2);
        end
    end
    
    valid_indices = (result_x(i, :) ~= 0) & (result_y(i, :) ~= 0);
    plot(result_x(i, valid_indices), result_y(i, valid_indices), ...
         'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'blue');
    
    if result_x(i, 1) ~= 0 && result_y(i, 1) ~= 0
        plot(result_x(i, 1), result_y(i, 1), 'go', 'MarkerSize', 6, ...
             'MarkerFaceColor', 'green', 'LineWidth', 2);
        text(result_x(i, 1), result_y(i, 1), '  龙头', 'FontSize', 8, ...
             'Color', 'green', 'FontWeight', 'bold');
    end
    
    last_valid = find(valid_indices, 1, 'last');
    if ~isempty(last_valid)
        plot(result_x(i, last_valid), result_y(i, last_valid), 'mo', ...
             'MarkerSize', 6, 'MarkerFaceColor', 'magenta', 'LineWidth', 2);
        text(result_x(i, last_valid), result_y(i, last_valid), '  龙尾', ...
             'FontSize', 8, 'Color', 'magenta', 'FontWeight', 'bold');
    end

    title(sprintf('时间: t = %.2f s ', ...
          t(i)), 'FontSize', 14, 'FontWeight', 'bold');

    all_x = [result_x(i, valid_indices), special_points_x];
    all_y = [result_y(i, valid_indices), special_points_y];
    if ~isempty(all_x)
        margin = 2;
        xlim([min(all_x) - margin, max(all_x) + margin]);
        ylim([min(all_y) - margin, max(all_y) + margin]);
    end
    
    hold off;
end
