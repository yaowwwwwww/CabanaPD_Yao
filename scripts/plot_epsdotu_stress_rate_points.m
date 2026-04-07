function plot_epsdotu_stress_rate_points()
  run_root = getenv('RUN_ROOT');
  if isempty(run_root)
    error('RUN_ROOT is not set');
  end

  data_dir = getenv('DATA_DIR');
  if isempty(data_dir)
    data_dir = fullfile(run_root, 'epsdotu_pointcloud_r5e-6');
  end

  tsv_path = fullfile(data_dir, 'stress_rate_points_by_epsdot_u.tsv');
  fid = fopen(tsv_path, 'r');
  if fid < 0
    error('Could not open %s', tsv_path);
  end
  cleaner = onCleanup(@() fclose(fid));
  textscan(fid, '%s', 1, 'Delimiter', '\n');
  C = textscan(fid, '%s%f%f%f%f', 'Delimiter', '\t');
  clear cleaner;

  labels = C{1};
  epsdot = C{4};
  ys = C{5};

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
  markers = {'o', 's', '^', 'x'};

  legend_labels = cell(numel(unique_labels), 1);
  for j = 1:numel(unique_labels)
    lbl = unique_labels{j};
    mask = strcmp(labels, lbl) & epsdot > 0.0 & ys > 0.0;
    x = epsdot(mask);
    y = ys(mask);
    color = colors(mod(j-1, size(colors,1)) + 1, :);
    marker = markers{mod(j-1, numel(markers)) + 1};
    plot(ax, x, y, marker, 'color', color, 'markersize', 5, ...
         'linewidth', 1.0, 'handlevisibility', 'off');

    x_fit = [];
    y_fit = [];
    if numel(x) >= 20
      x_min = min(x);
      x_max = max(x);
      if x_max > x_min
        edges = logspace(log10(x_min), log10(x_max), 28);
        for k = 1:(numel(edges) - 1)
          if k == numel(edges) - 1
            bin_mask = x >= edges(k) & x <= edges(k+1);
          else
            bin_mask = x >= edges(k) & x < edges(k+1);
          end
          if sum(bin_mask) >= 10
            x_fit(end+1,1) = median(x(bin_mask));
            y_fit(end+1,1) = median(y(bin_mask));
          end
        end
      end
    end

    if numel(x_fit) >= 2
      plot(ax, x_fit, y_fit, '-', 'color', color, 'linewidth', 3.0);
    end
    legend_labels{j} = sprintf('epsdot_u=%s 1/s', lbl);
  end

  grid(ax, 'on');
  box(ax, 'on');
  set(ax, 'fontsize', fs_tick, 'linewidth', 1.0, 'xscale', 'log');
  xlabel(ax, 'Plastic strain rate (1/s)', 'fontsize', fs_label);
  ylabel(ax, 'Flow stress (Pa)', 'fontsize', fs_label);
  title(ax, 'Near-impact point cloud: flow stress vs strain rate', 'fontsize', fs_title);
  legend(ax, legend_labels, 'location', 'eastoutside', 'fontsize', fs_legend);

  out_png = fullfile(data_dir, 'stress_rate_points_by_epsdot_u.png');
  print(fig, out_png, '-dpng', '-r220');
  close(fig);
  fprintf('saved %s\n', out_png);
end

plot_epsdotu_stress_rate_points();
