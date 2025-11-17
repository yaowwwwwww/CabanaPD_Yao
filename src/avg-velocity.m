files = dir("particles_*_all.csv");

% 按编号排序
[~, idx] = sort(cellfun(@(s) sscanf(s, "particles_%d_all.csv"), {files.name}));
files = files(idx);

outfile = fopen("avg_vmag.csv", "w");

for i = 1:length(files)
    fname = files(i).name;
    printf("Processing %s\n", fname);

    % 跳过表头
    data = csvread(fname, 1, 0);

    % 列索引（根据你的 CSV）
    col_type = 3;    % rank_0/type
    col_vx = 10;     % rank_0/velocities:0
    col_vy = 11;     % rank_0/velocities:1
    col_vz = 12;     % rank_0/velocities:2

    % 只要 type=1 的粒子
    mask = data(:, col_type) == 1;
    subset = data(mask, :);

    if isempty(subset)
        avg_vmag = NaN;
    else
        avg_vmag = mean(sqrt(sum(subset(:, [col_vx col_vy col_vz]).^2, 2)));
    end

    fprintf(outfile, "%.6e\n", avg_vmag);
end

fclose(outfile);
% 读取平均速度文件
data = load("avg_vmag.csv");

% 横坐标：时间步（0,1,2,...）
steps = (0:(length(data)-1))*0.5;

% 绘图
figure;
plot(steps, data, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
xlabel("Time(ns)");
ylabel("Average velocity magnitude(m/s)");
title("Evolution of average velocity (type=1 particles)");
 set(gca, "fontsize", 34);
grid on;
