function plot_epsdotu_stress_strain()
  run_root = getenv('RUN_ROOT');
  if isempty(run_root)
    error('RUN_ROOT is not set');
  end

  data_dir = getenv('DATA_DIR');
  if isempty(data_dir)
    data_dir = fullfile(run_root, 'epsdotu_tripanel');
  end

  tsv_path = fullfile(data_dir, 'mean_response_vs_time_by_epsdot_u.tsv');
  fid = fopen(tsv_path, 'r');
  if fid < 0
    error('Could not open %s', tsv_path);
  end
  cleaner = onCleanup(@() fclose(fid));
  textscan(fid, '%s', 1, 'Delimiter', '\n');
  C = textscan(fid, '%s%f%f%f%f%f%f%f', 'Delimiter', '\t');
  clear cleaner;

  labels = C{1};
  frame = C{2};
  mean_eps = C{4};
  mean_ys = C{6};

  if isempty(labels)
    error('No data rows found in %s', tsv_path);
  end

  fs_title = 22;
  fs_label = 20;
  fs_tick = 18;
  fs_legend = 16;

  fig = figure('visible', 'off');
  set(fig, 'position', [100 100 1200 900], 'color', 'w');
  ax = axes(fig);
  hold(ax, 'on');

  unique_labels = unique(labels, 'stable');
  colors = [
    0.80 0.18 0.18;
    0.10 0.45 0.80;
    0.92 0.58 0.08;
    0.25 0.65 0.35
  ];

  legend_labels = cell(numel(unique_labels), 1);
  for j = 1:numel(unique_labels)
    lbl = unique_labels{j};
    mask = strcmp(labels, lbl);
    [~, ord] = sort(frame(mask));
    epsx = mean_eps(mask); epsx = epsx(ord);
    ysx = mean_ys(mask); ysx = ysx(ord);
    color = colors(mod(j-1, size(colors,1)) + 1, :);
    plot(ax, epsx, ysx, 'linewidth', 2.2, 'color', color);
    legend_labels{j} = sprintf('epsdot_u=%s 1/s', lbl);
  end

  grid(ax, 'on');
  box(ax, 'on');
  set(ax, 'fontsize', fs_tick, 'linewidth', 1.0);
  xlabel(ax, 'Mean plastic strain', 'fontsize', fs_label);
  ylabel(ax, 'Mean yield stress (Pa)', 'fontsize', fs_label);
  title(ax, 'Mean yield stress vs mean plastic strain', 'fontsize', fs_title);
  legend(ax, legend_labels, 'location', 'eastoutside', 'fontsize', fs_legend);

  out_png = fullfile(data_dir, 'mean_stress_vs_mean_strain_by_epsdot_u.png');
  print(fig, out_png, '-dpng', '-r220');
  close(fig);
  fprintf('saved %s\n', out_png);
end

plot_epsdotu_stress_strain();
