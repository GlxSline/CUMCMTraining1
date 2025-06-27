%load data_2;

% result_v_max = max(result_v,[],2);

% figure;
% plot_x = linspace(1,100,1001);
% plot(plot_x, result_v_max);
% grid on;

result_v_max = max(result_v);

figure;
plot_x = linspace(1,5,5);
plot(plot_x, result_v_max);
grid on;


grid(max(result_v_max));
grid(2/max(result_v_max));

figure;
result_v4 = result_v(:,4);
plot_x = linspace(1,100,10001);
plot(plot_x, result_v4);
grid on;