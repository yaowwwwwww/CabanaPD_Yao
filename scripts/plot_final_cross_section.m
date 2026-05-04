% Plot the final-frame x-z cross-section through the impact center.
%
% Usage:
%   CASE_DIR=/path/to/case octave-cli --quiet scripts/plot_final_cross_section.m
%
% Optional environment variables:
%   SLICE_WIDTH_UM  full slice thickness in y, default 1.0 um
%   OUT_PNG         output image path
%   OUT_CSV         output cross-section CSV path

case_dir = getenv("CASE_DIR");
if isempty(case_dir)
    case_dir = pwd;
end

slice_width_um_env = getenv("SLICE_WIDTH_UM");
if isempty(slice_width_um_env)
    slice_width_um = 1.0;
else
    slice_width_um = str2double(slice_width_um_env);
end
slice_half_width_m = 0.5 * slice_width_um * 1.0e-6;

out_png = getenv("OUT_PNG");
if isempty(out_png)
    out_png = fullfile(case_dir, "final_cross_section.png");
end

out_csv = getenv("OUT_CSV");
if isempty(out_csv)
    out_csv = fullfile(case_dir, "final_cross_section_slice.csv");
end

files = dir(fullfile(case_dir, "particles_*_all.csv"));
if isempty(files)
    error("No particles_*_all.csv files found in %s", case_dir);
end

frame_ids = zeros(numel(files), 1);
for i = 1:numel(files)
    tok = regexp(files(i).name, "particles_(\\d+)_all\\.csv", "tokens");
    if isempty(tok)
        frame_ids(i) = -1;
    else
        frame_ids(i) = str2double(tok{1}{1});
    end
end
[~, order] = sort(frame_ids);
files = files(order);
frame_ids = frame_ids(order);

final_file = fullfile(case_dir, files(end).name);
final_frame = frame_ids(end);

fid = fopen(final_file, "r");
if fid < 0
    error("Cannot open %s", final_file);
end
header_line = fgetl(fid);
fclose(fid);
headers = strsplit(strrep(header_line, "\"", ""), ",");

col_type = find(strcmp(headers, "rank_0/type"), 1);
col_x = find(strcmp(headers, "Points:0"), 1);
col_y = find(strcmp(headers, "Points:1"), 1);
col_z = find(strcmp(headers, "Points:2"), 1);

if isempty(col_type) || isempty(col_x) || isempty(col_y) || isempty(col_z)
    error("Missing required columns in %s", final_file);
end

data = csvread(final_file, 1, 0);

ball_type = 0;
substrate_type = 1;
ball = data(data(:, col_type) == ball_type, :);
substrate = data(data(:, col_type) == substrate_type, :);

if isempty(ball)
    warning("No projectile/type=0 points found; using y=0 for slice center.");
    y_center = 0.0;
else
    y_center = median(ball(:, col_y));
end

slice_mask = abs(data(:, col_y) - y_center) <= slice_half_width_m;
slice_data = data(slice_mask, :);
if isempty(slice_data)
    error("Slice is empty. Increase SLICE_WIDTH_UM. Current width = %.6g um", slice_width_um);
end

slice_ball = slice_data(slice_data(:, col_type) == ball_type, :);
slice_sub = slice_data(slice_data(:, col_type) == substrate_type, :);

fid = fopen(out_csv, "w");
if fid < 0
    error("Cannot write %s", out_csv);
end
fprintf(fid, "type,x_um,y_um,z_um\n");
for i = 1:size(slice_data, 1)
    fprintf(fid, "%.0f,%.9g,%.9g,%.9g\n", slice_data(i, col_type), ...
            slice_data(i, col_x) * 1.0e6, slice_data(i, col_y) * 1.0e6, ...
            slice_data(i, col_z) * 1.0e6);
end
fclose(fid);

gp_file = [tempname(), ".gp"];
fid = fopen(gp_file, "w");
if fid < 0
    error("Cannot write temporary gnuplot script");
end
fprintf(fid, "set terminal pngcairo size 1500,1000 enhanced font 'Arial,24'\n");
fprintf(fid, "set output '%s'\n", out_png);
fprintf(fid, "set datafile separator comma\n");
fprintf(fid, "set key outside right top box opaque\n");
fprintf(fid, "set grid\n");
fprintf(fid, "set border linewidth 1.5\n");
fprintf(fid, "set size ratio -1\n");
fprintf(fid, "set xlabel 'Distance x ({/Symbol m}m)'\n");
fprintf(fid, "set ylabel 'Height z ({/Symbol m}m)'\n");
fprintf(fid, "set title 'Final frame %d x-z cross-section (y into page)'\n", final_frame);
fprintf(fid, "plot '%s' using (($1==1)?$2:1/0):4 every ::1 with points pt 7 ps 0.55 lc rgb '#1f5fd0' title 'Substrate', \\\n", out_csv);
fprintf(fid, "     '%s' using (($1==0)?$2:1/0):4 every ::1 with points pt 7 ps 0.75 lc rgb '#e31a1c' title 'Projectile', \\\n", out_csv);
fprintf(fid, "     0 with lines dt 2 lw 2 lc rgb 'black' title 'Initial surface z=0'\n");
fclose(fid);

[status, output] = system(sprintf("gnuplot '%s'", gp_file));
unlink(gp_file);
if status ~= 0
    error("gnuplot failed: %s", output);
end

fprintf("Final frame: %s\n", final_file);
fprintf("Slice center y: %.9g um\n", y_center * 1.0e6);
fprintf("Slice width: %.9g um\n", slice_width_um);
fprintf("Projectile points in slice: %d\n", size(slice_ball, 1));
fprintf("Substrate points in slice: %d\n", size(slice_sub, 1));
fprintf("Wrote image: %s\n", out_png);
fprintf("Wrote slice CSV: %s\n", out_csv);
