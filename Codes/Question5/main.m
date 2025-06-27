load data;

result_v_max = max(result_v,[],2);

figure;
plot((0:100),result_v_max);
grid on;

