% recomputehmax.m
% Output:
% case_name, vin, vout, CoR, Lateralmax,
% h_surface(m), h_penetration(m), A_residual(m2)

clear; clc;

root_dir = pwd;

summary_env = getenv("SUMMARY_FILE");
if isempty(summary_env)
    logfile = fullfile(root_dir, "summary_depth_recomputed.txt");
else
    logfile = summary_env;
end

single_case_dir = getenv("SINGLE_CASE_DIR");
single_case = ~isempty(single_case_dir);

% ---------- open summary file ----------
if exist(logfile, "file") == 2
    fid_log = fopen(logfile, "a");   % append
else
    fid_log = fopen(logfile, "w");   % create new file
    fprintf(fid_log, "# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_surface(m)    h_penetration(m)    A_residual(m2)\n");
end
 
if fid_log < 0
    error("Cannot create logfile: %s", logfile);
end
 
if single_case
    [case_path, base, ext] = fileparts(single_case_dir);
    case_name = [base ext];
    dirs = struct("name", case_name, "folder", case_path);
else
    dirs = dir("v_*");
    dirs = dirs([dirs.isdir]);

    if isempty(dirs)
        dirs = dir("run_vin_*");
        dirs = dirs([dirs.isdir]);
    end
end

if isempty(dirs)
    error("No v_* or run_vin_* directories found.");
end

fprintf("Found %d cases.\n", numel(dirs));

for id = 1:numel(dirs)

    case_name = dirs(id).name;
    fprintf("\n==============================\n");
    fprintf("Case %d / %d : %s\n", id, numel(dirs), case_name);

    if single_case
        case_dir = single_case_dir;
    else
        case_dir = fullfile(root_dir, case_name);
    end

    cd(case_dir);

    files = dir("particles_*_all.csv");
    if isempty(files)
        fprintf("  [WARN] No particles_*_all.csv, skip.\n");
        fprintf(fid_log, "%s NaN NaN NaN NaN NaN NaN NaN\n", case_name);
        cd(root_dir);
        continue;
    end

    [~, idx] = sort(cellfun(@(s) sscanf(s, "particles_%d_all.csv"), {files.name}));
    files = files(idx);
    nFrames = numel(files);

    fid_hdr = fopen(files(1).name, "r");
    header_line = fgetl(fid_hdr);
    fclose(fid_hdr);

    headers = strsplit(strrep(header_line, "\"", ""), ",");

    col_type = find(strcmp(headers, "rank_0/type"), 1);
    col_vx   = find(strcmp(headers, "rank_0/velocities:0"), 1);
    col_vy   = find(strcmp(headers, "rank_0/velocities:1"), 1);
    col_vz   = find(strcmp(headers, "rank_0/velocities:2"), 1);
    col_x    = find(strcmp(headers, "Points:0"), 1);
    col_y    = find(strcmp(headers, "Points:1"), 1);
    col_z    = find(strcmp(headers, "Points:2"), 1);

    if isempty(col_type) || isempty(col_vx) || isempty(col_vy) || isempty(col_vz) || isempty(col_x) || isempty(col_y) || isempty(col_z)
        col_type = 3;
        col_vx = 10; col_vy = 11; col_vz = 12;
        col_x = 13;  col_y = 14;  col_z = 15;
    end

    ball_type = 0;
    substrate_type = 1;

    % ==========================================================
    % vin / vout / CoR / Lateralmax
    % ==========================================================
    data_first = csvread(files(1).name, 1, 0);
    first_projectile = data_first(data_first(:, col_type) == ball_type, :);

    if isempty(first_projectile)
        d0 = NaN;
    else
        dx0 = max(first_projectile(:, col_x)) - min(first_projectile(:, col_x));
        dy0 = max(first_projectile(:, col_y)) - min(first_projectile(:, col_y));
        d0 = 0.5 * (dx0 + dy0);
    end

    vavg  = nan(nFrames,1);
    dcoef = nan(nFrames,1);

    for k = 1:nFrames
        data = csvread(files(k).name, 1, 0);
        projectile = data(data(:, col_type) == ball_type, :);

        if isempty(projectile)
            continue;
        end

        vavg(k) = mean(projectile(:, col_vz));

        dx = max(projectile(:, col_x)) - min(projectile(:, col_x));
        dy = max(projectile(:, col_y)) - min(projectile(:, col_y));
        dz = abs(max(projectile(:, col_z)) - min(projectile(:, col_z)));

        lateral = 0.5 * (dx + dy);

        if isnan(d0) || d0 <= 0
            dcoef(k) = NaN;
        else
            dcoef(k) = (lateral - dz) / d0;
        end
    end

    vin = vavg(1);
    vout = vavg(end);
    CoR = vout / vin;
    Lateralmax = max(dcoef);

    fprintf("  vin = %.4f, vout = %.4f, CoR = %.4f, Lateralmax = %.4f\n", ...
            vin, vout, CoR, Lateralmax);

    % ==========================================================
    % Depth calculations
    % ==========================================================
    data0 = csvread(files(1).name, 1, 0);
    sub0 = data0(data0(:, col_type) == substrate_type, :);

    if isempty(sub0)
        all_types = unique(data0(:, col_type));
        cand = all_types(all_types ~= ball_type);
        if ~isempty(cand)
            substrate_type = cand(1);
            sub0 = data0(data0(:, col_type) == substrate_type, :);
        end
    end

    if isempty(sub0)
        warning("No substrate particles found.");
        h_surface = NaN;
        h_penetration = NaN;
        A_residual = NaN;
    else
        BINS = 300;
        TOPK = 15;
        DEPTH_K = 20;
        BAND_FRAC = 0.01;

        y_mid = median(sub0(:, col_y));
        eps_band = (max(sub0(:, col_y)) - min(sub0(:, col_y))) * BAND_FRAC;

        line0 = sub0(abs(sub0(:, col_y) - y_mid) <= eps_band, :);
        if isempty(line0)
            line0 = sub0;
        end

        xmin = min(sub0(:, col_x));
        xmax = max(sub0(:, col_x));
        edges = linspace(xmin, xmax, BINS + 1);
        dx_bin = (xmax - xmin) / BINS;

        % ---------- initial surface envelope ----------
        z0_env = nan(1, BINS);

        for b = 1:BINS
            xl = edges(b);
            xr = edges(b+1);

            m0 = (line0(:, col_x) >= xl) & (line0(:, col_x) < xr);

            if any(m0)
                zlist = sort(line0(m0, col_z), "descend");
                kk = min(TOPK, numel(zlist));
                z0_env(b) = mean(zlist(1:kk));
            end
        end

        % ---------- final frame ----------
        dataN = csvread(files(end).name, 1, 0);
        subN = dataN(dataN(:, col_type) == substrate_type, :);

        if isempty(subN)
            h_surface = NaN;
            h_penetration = NaN;
            A_residual = NaN;
        else
            lineN = subN(abs(subN(:, col_y) - y_mid) <= eps_band, :);
            if isempty(lineN)
                lineN = subN;
            end

            % ---------- final surface envelope ----------
            zN_env = nan(1, BINS);

            for b = 1:BINS
                xl = edges(b);
                xr = edges(b+1);

                mN = (lineN(:, col_x) >= xl) & (lineN(:, col_x) < xr);

                if any(mN)
                    zlist = sort(lineN(mN, col_z), "descend");
                    kk = min(TOPK, numel(zlist));
                    zN_env(b) = mean(zlist(1:kk));
                end
            end

            valid = ~isnan(z0_env) & ~isnan(zN_env);

            if any(valid)
                h_full = zeros(1, BINS);
                h_tmp = z0_env(valid) - zN_env(valid);
                h_tmp(h_tmp < 0) = 0;
                h_full(valid) = h_tmp;

                h_vals = h_full(h_full > 0);

                if isempty(h_vals)
                    h_surface = 0;
                else
                    h_sorted = sort(h_vals, "descend");
                    K = min(DEPTH_K, numel(h_sorted));
                    h_surface = mean(h_sorted(1:K));
                end

                A_residual = sum(h_full) * dx_bin;
            else
                h_surface = NaN;
                A_residual = NaN;
            end

 
% Penetration depth from free-surface envelope
% ======================================================

h_vals = h_full(h_full > 0);

if isempty(h_vals)
    h_penetration = 0;
else
    h_penetration = max(h_vals);
end
        end
    end

    fprintf("  h_surface      = %.6e m\n", h_surface);
    fprintf("  h_penetration  = %.6e m\n", h_penetration);
    fprintf("  A_residual     = %.6e m^2\n", A_residual);

    summary_name = case_name;
    if length(summary_name) >= 4 && strcmp(summary_name(1:4), "run_")
        summary_name = summary_name(5:end);
    end

    fprintf(fid_log, "%s %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n", ...
            summary_name, vin, vout, CoR, Lateralmax, ...
            h_surface, h_penetration, A_residual);

    fflush(fid_log);
    cd(root_dir);
end

fclose(fid_log);
fprintf("\nDone. Summary saved to %s\n", logfile);