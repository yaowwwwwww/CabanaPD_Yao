% plot_compare_lit_vs_mine.m
clear; clc;

% ---------- input ----------
csvfile = "983ada83-fa85-4f9e-95e4-d8f7196721bc.csv";  % <-- change if needed

fid = fopen(csvfile, "r");
if fid < 0
    error("Cannot open file: %s", csvfile);
end

% ---------- read header ----------
header = fgetl(fid); %#ok<NASGU>

% ---------- storage ----------
% Literature (already in GPa)
sr_lit = [];
H_lit  = [];

% Yours: ys (Pa), strain rate (1/s), hardness (Pa)
ys_all = [];
sr_my  = [];
Hd_my  = [];

% ---------- parse lines ----------
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    line = strtrim(line);
    if isempty(line), continue; end

    parts = strsplit(line, ",");

    % need at least 4 comma-separated fields
    if numel(parts) < 4
        continue;
    end

    % --- literature cols ---
    sr1 = str2double(strtrim(parts{1}));
    H1  = str2double(strtrim(parts{2}));

    if isfinite(sr1) && isfinite(H1)
        sr_lit(end+1,1) = sr1;
        H_lit(end+1,1)  = H1;   % GPa
    end

    % --- your packed col: "ys  strain_rate  h_d" ---
    packed = strtrim(parts{4});
    vals = sscanf(packed, "%e %e %e");  % [ys, sr, h_d]
    if numel(vals) == 3 && all(isfinite(vals))
        ys_all(end+1,1) = vals(1);   % Pa
        sr_my(end+1,1)  = vals(2);   % 1/s
        Hd_my(end+1,1)  = vals(3);   % Pa
    end
end

fclose(fid);

if isempty(sr_lit) || isempty(sr_my)
    error("Parsed data seems empty. Check file formatting.");
end

% Convert your hardness to GPa for comparison
H_my = Hd_my / 1e9;

% ---------- plot ----------
figure(1); clf; hold on;

LW = 2.0;
MS = 9;
FS = 28;

% Literature scatter (gray circles)
plot(sr_lit, H_lit, 'o', ...
     'LineStyle','none', 'MarkerSize', MS, 'LineWidth', LW);

% Your data grouped by ys (scatter only)
ys_set = unique(ys_all);
line_styles = {'rs','b^','gd','m+','cx','kp','y*','ro','bs'}; % markers only

h_list = {};
legend_str = {};

% First legend entry: literature
h_list{end+1} = plot(NaN, NaN, 'o', 'LineStyle','none', 'MarkerSize', MS, 'LineWidth', LW);
legend_str{end+1} = "Literature";

for k = 1:numel(ys_set)
    ys_k = ys_set(k);
    idx  = (ys_all == ys_k);

    sr_k = sr_my(idx);
    H_k  = H_my(idx);

    good = isfinite(sr_k) & isfinite(H_k) & (sr_k > 0) & (H_k > 0);
    sr_k = sr_k(good);
    H_k  = H_k(good);

    style = line_styles{mod(k-1, numel(line_styles)) + 1};

    plot(sr_k, H_k, style, ...
         'LineStyle','none', 'MarkerSize', MS, 'LineWidth', LW);

    % Add legend handle (dummy point) so legend shows one per ys
    h_list{end+1} = plot(NaN, NaN, style, 'LineStyle','none', 'MarkerSize', MS, 'LineWidth', LW);
    legend_str{end+1} = sprintf('My data: ys_{Cu}=%.2g GPa', ys_k/1e9);
end

set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');

xlabel('Strain rate (s^{-1})');
ylabel('Hardness (GPa)');
title('Hardness vs Strain rate: Literature vs My data');
legend([h_list{:}], legend_str, 'Location', 'northwest');

grid on;
set(gca, 'fontsize', FS);

% Optional: set x/y limits if you want
% xlim([1e4 1e8]);
% ylim([0.5 20]);

fprintf("Parsed %d literature points, %d my points.\n", numel(sr_lit), numel(sr_my));

