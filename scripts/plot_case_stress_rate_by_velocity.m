manifest = getenv('CASE_HISTORY_MANIFEST');
if isempty(manifest)
  error('CASE_HISTORY_MANIFEST is not set');
end

fid = fopen(manifest, 'r');
if fid < 0
  error('cannot open manifest: %s', manifest);
end
fgetl(fid);
M = textscan(fid, '%s %s %s %s %f %f %f', 'Delimiter', '\t');
fclose(fid);

case_names = M{1};
raw_tsvs = M{2};

vel_tags = {};
for i = 1:numel(case_names)
  tok = regexp(case_names{i}, '^(v\d+)_', 'tokens');
  if ~isempty(tok)
    vel_tags{end + 1} = tok{1}{1}; %#ok<AGROW>
  end
end
vel_tags = unique(vel_tags);

for v = 1:numel(vel_tags)
  vel = vel_tags{v};
  fig = figure('visible', 'off');
  hold on;
  plotted = 0;

  for i = 1:numel(case_names)
    if isempty(regexp(case_names{i}, ['^' vel '_'], 'once'))
      continue;
    end
    tok = regexp(case_names{i}, '_K0_([^_]+)_', 'tokens');
    if isempty(tok)
      continue;
    end
    k0 = tok{1}{1};
    raw = dlmread(raw_tsvs{i}, '\t', 1, 0);
    if rows(raw) < 6
      continue;
    end
    x = raw(:, 2);
    y = raw(:, 3);
    valid = (x > 0) & (y > 0);
    x = x(valid);
    y = y(valid);
    if numel(x) < 6
      continue;
    end

    logx = log10(x);
    nbins = min(30, max(12, floor(numel(x) / 20)));
    edges = linspace(min(logx), max(logx), nbins + 1);
    bx = [];
    by = [];
    for b = 1:nbins
      if b < nbins
        idx = find(logx >= edges(b) & logx < edges(b + 1));
      else
        idx = find(logx >= edges(b) & logx <= edges(b + 1));
      end
      if numel(idx) < 3
        continue;
      end
      bx(end + 1, 1) = 10 ^ mean(logx(idx)); %#ok<AGROW>
      by(end + 1, 1) = mean(y(idx)); %#ok<AGROW>
    end
    if numel(bx) < 2
      continue;
    end
    [bx, ord] = sort(bx);
    by = by(ord);
    semilogx(bx, by ./ 1e9, '-o', 'LineWidth', 1.8, 'MarkerSize', 4.5, ...
             'MarkerFaceColor', 'none', 'DisplayName', ['K0=' k0]);
    plotted = plotted + 1;
  end

  if plotted == 0
    close(fig);
    continue;
  end

  xlabel('Plastic strain rate (1/s)');
  ylabel('Flow stress (GPa)');
  title([vel ': mean flow stress by strain-rate bin']);
  grid on;
  box on;
  legend('show', 'Location', 'northwest');

  [out_dir, ~, ~] = fileparts(manifest);
  out_png = fullfile(out_dir, [vel '_stress_rate_by_k0.png']);
  print(fig, out_png, '-dpng', '-r220');
  close(fig);
end
