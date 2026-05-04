% recomputehmax.m
% 在当前目录下，遍历所有 run_vin_* 子目录，
% 对每个 case 重新计算：vin, vout, CoR, Lateralmax, h_max,
% h_residual, A_residual, h_mean_residual，并写一个总 summary 文件。

clear; clc;

root_dir = pwd;

% ---------- 输出汇总日志 ----------
summary_env = getenv("SUMMARY_FILE");
if isempty(summary_env)
    logfile = fullfile(root_dir, "summary_cor_hmax_recomputed.txt");
else
    logfile = summary_env;
end

single_case_dir = getenv("SINGLE_CASE_DIR");
single_case = ~isempty(single_case_dir);

if single_case
    % append-only mode for a single case (called by scan script)
    if exist(logfile, "file") == 2
        fid_log = fopen(logfile, "a");
    else
        fid_log = fopen(logfile, "w");
        if fid_log >= 0
            fprintf(fid_log, "# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_max(m)    best_frame    h_residual(m)    A_residual(m2)    h_mean_residual(m)\n");
        end
    end
else
    fid_log = fopen(logfile, "w");
    if fid_log >= 0
        fprintf(fid_log, "# case_name    vin(m/s)    vout(m/s)    CoR    Lateralmax    h_max(m)    best_frame    h_residual(m)    A_residual(m2)    h_mean_residual(m)\n");
    end
end

if fid_log < 0
    error("无法创建日志文件: %s", logfile);
end

% ---------- 找所有模拟子目录 ----------
if single_case
    % compute only the provided case; preserve full name with decimals
    [case_path, base, ext] = fileparts(single_case_dir);
    case_name = [base ext];
    dirs = struct("name", case_name, "folder", case_path);
else
    dirs = dir("run_vin_*");
    dirs = dirs([dirs.isdir]);
end

if isempty(dirs)
    error("当前目录下没有找到 run_vin_* 子目录");
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

    % ----------------------------------------------------------
    % 读取所有粒子 CSV 帧
    % ----------------------------------------------------------
    files = dir("particles_*_all.csv");
    if isempty(files)
        fprintf("  [WARN] No particles_*_all.csv in %s, skip.\n", case_name);
        fprintf(fid_log, "%s  NaN  NaN  NaN  NaN  NaN  -1  NaN  NaN  NaN\n", case_name);
        fflush(fid_log);
        cd(root_dir);
        continue;
    end

    % 按帧号排序
    [~, idx] = sort( cellfun(@(s) sscanf(s, "particles_%d_all.csv"), {files.name}) );
    files = files(idx);
    nFrames = numel(files);

    % --------- 列索引（从 CSV 表头自动匹配） ---------
    fid_hdr = fopen(files(1).name, "r");
    if fid_hdr < 0
        error("Cannot open CSV header: %s", files(1).name);
    end
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
        % Fallback: no/unknown header, use legacy fixed columns
        col_type = 3;
        col_vx = 10; col_vy = 11; col_vz = 12;
        col_x = 13;  col_y = 14;  col_z = 15;
    end

% ==========================================================
% 1) 计算 vin / vout / CoR / Lateralmax （用你原来的公式）
%    projectile: type = ball_type
% ==========================================================
ball_type = 0;
substrate_type = 1;
    vavg  = nan(nFrames,1);
    vmag  = nan(nFrames,1);
    dcoef = nan(nFrames,1);

    for k = 1:nFrames
        fname = files(k).name;
        data = csvread(fname, 1, 0);

        % select projectile particles
        mask = data(:, col_type) == ball_type;
        subset = data(mask, :);

        if isempty(subset)
            avg_vz = NaN;
            avg_vmag = NaN;
            deform_coeff = NaN;
        else
            % average normal velocity (signed)
            avg_vz = mean(subset(:, col_vz));
            avg_vmag = mean(sqrt(sum(subset(:, [col_vx col_vy col_vz]).^2, 2)));

            % ---- lateral deformation coefficient ----
            dx = max(subset(:, col_x)) - min(subset(:, col_x));
            dy = max(subset(:, col_y)) - min(subset(:, col_y));
            dz = abs(max(subset(:, col_z)) - min(subset(:, col_z)));
            d0 = 12e-6;               % initial diameter
            lateral = (dx + dy)/2;    % mean lateral span
            deform_coeff = (lateral - dz) / d0;
        end

        vavg(k)  = avg_vz;
        vmag(k)  = avg_vmag;
        dcoef(k) = deform_coeff;
    end

    % ---- record avg_vmag per frame (same format as old avg-velocity.m) ----
    fout_vmag = fopen("avg_vmag.csv", "w");
    if fout_vmag >= 0
        for k = 1:nFrames
            fprintf(fout_vmag, "%.6e\n", vmag(k));
        end
        fclose(fout_vmag);
    end

    % --------- simplest CoR: last frame / first frame ---------
    if all(isnan(vavg))
        vin  = NaN;
        vout = NaN;
        CoR  = NaN;
    else
        vin  = vavg(1);          % 初始速度（第一帧）
        vout = vavg(end);        % 最后一帧速度

        CoR = vout / vin;   % 恢复系数
    end

    Lateralmax = max(dcoef);

    fprintf("  vin = %.4f m/s, vout = %.4f m/s, CoR = %.4f, Lateralmax = %.4f\n", ...
            vin, vout, CoR, Lateralmax);

    % ==========================================================
    % 2) h_max + 残余压痕（substrate = substrate_type）
    % ==========================================================

    % --- baseline surface from frame 1 (substrate only)
    data0 = csvread(files(1).name, 1, 0);
    sub0  = data0(data0(:,col_type)==substrate_type, :);

    % 如果 type=1 根本不存在，自动猜一个“非 projectile”的类型做基底
    if isempty(sub0)
        all_types = unique(data0(:,col_type));
        cand = all_types(all_types ~= ball_type); % 排除 projectile 类型
        if isempty(cand)
            warning('No substrate points found in frame 1. h_max analysis skipped.');
            h_max = NaN;
            h_idx = -1;
            h_residual      = NaN;
            A_residual      = NaN;
            h_mean_residual = NaN;

            fprintf(fid_log, "%s  %.6e  %.6e  %.6e  %.6e  %.6e  %d  %.6e  %.6e  %.6e\n", ...
                    case_name, vin, vout, CoR, Lateralmax, h_max, h_idx, ...
                    h_residual, A_residual, h_mean_residual);
            fflush(fid_log);
            cd(root_dir);
            continue;
        else
            substrate_type = cand(1);
            sub0  = data0(data0(:,col_type)==substrate_type, :);
        end
    end

    if isempty(sub0)
        warning('Still no substrate points after auto-detect. h_max NaN.');
        h_max = NaN;
        h_idx = -1;
        h_residual      = NaN;
        A_residual      = NaN;
        h_mean_residual = NaN;
    else
        % parameters for mid-line Top-K (still robust, but only for one frame)
        BINS      = 300;
        TOPK      = 15;
        MINPTS    = 3;
        BAND_FRAC = 0.01;

        % 中线带
        y_mid    = median(sub0(:,col_y));
        eps_band = (max(sub0(:,col_y)) - min(sub0(:,col_y))) * BAND_FRAC;
        line0    = sub0(abs(sub0(:,col_y) - y_mid) <= eps_band, :);
        if isempty(line0)
            line0 = sub0;
        end

        % x 区间与 baseline envelope
        xmin  = min(sub0(:,col_x));
        xmax  = max(sub0(:,col_x));
        edges = linspace(xmin, xmax, BINS+1);

        z0_env = nan(1,BINS);
        for b = 1:BINS
            xl = edges(b); xr = edges(b+1);
            m0 = (line0(:,col_x) >= xl) & (line0(:,col_x) < xr);
            if any(m0)
                zlist = sort(line0(m0, col_z), 'descend');
                kk = min(max(MINPTS, TOPK), numel(zlist));
                z0_env(b) = mean(zlist(1:kk));
            end
        end

        % --------- use the final frame for the reported indentation depth ---------
        h_idx = numel(files);

        % --------- 在这一帧上计算压痕 ---------
        data_k = csvread(files(h_idx).name, 1, 0);
        sub_k  = data_k(data_k(:,col_type)==substrate_type, :);
        if isempty(sub_k)
            warning('No substrate points in selected frame. h_max NaN.');
            h_max = NaN;
            h_idx = -1;
        else
            lineN  = sub_k(abs(sub_k(:,col_y) - y_mid) <= eps_band, :);
            if isempty(lineN)
                lineN = sub_k;
            end

            zN_env = nan(1,BINS);
            for b = 1:BINS
                xl = edges(b); xr = edges(b+1);
                mN = (lineN(:,col_x) >= xl) & (lineN(:,col_x) < xr);
                if any(mN)
                    zlist = sort(lineN(mN, col_z), 'descend');
                    kk = min(max(MINPTS, TOPK), numel(zlist));
                    zN_env(b) = mean(zlist(1:kk));
                end
            end

            valid = ~isnan(z0_env) & ~isnan(zN_env);

            if ~any(valid)
                % 如果 Top-K 包络完全失败，退回用全局“最高表面-最低表面”作为近似
                z0_top = max(sub0(:,col_z));
                zN_min = min(sub_k(:,col_z));
                h_max  = max(0, z0_top - zN_min);
            else
                h_k = z0_env(valid) - zN_env(valid);
                h_k(h_k < 0) = 0;     % 防止数值噪声导致“负压痕”

                h_vals = h_k(h_k > 0);
                if isempty(h_vals)
                    h_max = 0;
                else
                    h_sorted = sort(h_vals, 'descend');
                    K = min(20, numel(h_sorted));
                    h_max = mean(h_sorted(1:K));
                end
            end
        end
    end

    fprintf('  h_max = %.6e (frame index = %d)\n', h_max, h_idx);

    % ======================================================================
    % 残余压痕 h_residual / A_residual / h_mean_residual 计算
    % 使用同一个 baseline (z0_env / xmin / xmax / y_mid / eps_band)
    % ======================================================================
    if ~exist('z0_env', 'var') || all(isnan(z0_env))
        warning('z0_env all NaN, residual indentation set to NaN.');
        h_residual      = NaN;
        A_residual      = NaN;
        h_mean_residual = NaN;
    else
        nFrames      = numel(files);
        h_series     = nan(nFrames, 1);  % 瞬时特征深度（Top-K）
        A_series     = nan(nFrames, 1);  % 截面面积
        hmean_series = nan(nFrames, 1);  % 等效平均深度

        dx_bin = (max(sub0(:,col_x)) - min(sub0(:,col_x))) / BINS;

        for k = 1:nFrames
            data_k = csvread(files(k).name, 1, 0);
            sub_k  = data_k(data_k(:,col_type)==substrate_type, :);
            if isempty(sub_k)
                h_series(k)     = NaN;
                A_series(k)     = NaN;
                hmean_series(k) = NaN;
                continue;
            end

            % 当前帧中线带
            lineN = sub_k(abs(sub_k(:,col_y) - y_mid) <= eps_band, :);
            if isempty(lineN)
                lineN = sub_k;
            end

            % 当前帧 envelope zN_env
            zN_env = nan(1,BINS);
            for b = 1:BINS
                xl = edges(b); xr = edges(b+1);
                mN = (lineN(:,col_x) >= xl) & (lineN(:,col_x) < xr);
                if any(mN)
                    zlist = sort(lineN(mN, col_z), 'descend');
                    kk = min(max(MINPTS, TOPK), numel(zlist));
                    zN_env(b) = mean(zlist(1:kk));
                end
            end

            valid = ~isnan(z0_env) & ~isnan(zN_env);

            if ~any(valid)
                % 兜底：最高初始表面 - 最低当前表面
                z0_top = max(sub0(:,col_z));
                zN_min = min(sub_k(:,col_z));
                h_frame = max(0, z0_top - zN_min);

                A_frame      = h_frame * (max(sub0(:,col_x)) - min(sub0(:,col_x)));
                h_mean_frame = h_frame;
            else
                % 完整 h(x) 序列（无效 bin 视为 0）
                h_full = zeros(1, BINS);
                h_k    = z0_env(valid) - zN_env(valid);
                h_k(h_k < 0) = 0;
                h_full(valid) = h_k;

                % 1) 瞬时特征深度：Top-K 的平均
                h_vals = h_full(h_full > 0);
                if isempty(h_vals)
                    h_frame = 0;
                else
                    h_sorted = sort(h_vals, 'descend');
                    K = min(20, numel(h_sorted));
                    h_frame = mean(h_sorted(1:K));
                end

                % 2) 截面面积：∫ h(x) dx
                A_frame = sum(h_full) * dx_bin;

                % 3) 等效平均深度：面积 / 有效坑宽
                if any(h_full > 0)
                    h_thr     = 0.1 * max(h_full);      % 阈值：10% 最大深度
                    mask      = h_full > h_thr;
                    L_crater  = sum(mask) * dx_bin;
                    if L_crater > 0
                        h_mean_frame = A_frame / L_crater;
                    else
                        h_mean_frame = 0;
                    end
                else
                    h_mean_frame = 0;
                end
            end

            h_series(k)     = h_frame;
            A_series(k)     = A_frame;
            hmean_series(k) = h_mean_frame;
        end

        % ---- residual indentation: use final frame only ----
        nFrames_valid = sum(~isnan(h_series));
        if nFrames_valid == 0
            h_residual      = NaN;
            A_residual      = NaN;
            h_mean_residual = NaN;
        else
            h_residual      = h_series(end);
            A_residual      = A_series(end);
            h_mean_residual = hmean_series(end);
        end
    end

    fprintf("Residual indentation metrics (last frame):\n");
    fprintf("  h_residual      = %.6e m\n",  h_residual);
    fprintf("  A_residual      = %.6e m^2\n", A_residual);
    fprintf("  h_mean_residual = %.6e m\n",  h_mean_residual);

    % ---------- 写入一行 summary ----------
    summary_name = case_name;
    if length(summary_name) >= 4 && strcmp(summary_name(1:4), "run_")
        summary_name = summary_name(5:end);
    end
    fprintf(fid_log, "%s  %.6e  %.6e  %.6e  %.6e  %.6e  %d  %.6e  %.6e  %.6e\n", ...
            summary_name, vin, vout, CoR, Lateralmax, h_max, h_idx, ...
            h_residual, A_residual, h_mean_residual);
    fflush(fid_log);

    % 回到根目录，继续下一个 case
    cd(root_dir);
end

fclose(fid_log);
fprintf("\nAll done. Summary saved to %s\n", logfile);
