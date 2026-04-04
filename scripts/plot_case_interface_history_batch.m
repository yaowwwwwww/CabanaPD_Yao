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
mean_tsvs = M{3};
out_pngs = M{4};

for i = 1:numel(case_names)
  raw_tsv = raw_tsvs{i};
  out_png = out_pngs{i};

  raw = dlmread(raw_tsv, '\t', 1, 0);
  fig = figure('visible', 'off');
  hold on;
  h_scatter = [];
  h_curve = [];

  if rows(raw) > 0
    stride = max(1, floor(rows(raw) / 4000));
    idx = 1:stride:rows(raw);
    h_scatter = plot(raw(idx, 2), raw(idx, 3) ./ 1e9, '.', ...
                     'Color', [0.6 0.6 0.6], 'MarkerSize', 4);
  end

  % Red curve: mean flow stress in strain-rate bins.
  if rows(raw) > 5
    mx = raw(:, 2);
    my = raw(:, 3);
    validm = (mx > 0) & (my > 0);
    mx = mx(validm);
    my = my(validm);
    if numel(mx) > 5
      logx = log10(mx);
      nbins = min(30, max(12, floor(numel(mx) / 20)));
      edges = linspace(min(logx), max(logx), nbins + 1);
      bx = [];
      by = [];
      for b = 1:nbins
        if b < nbins
          idxb = find(logx >= edges(b) & logx < edges(b + 1));
        else
          idxb = find(logx >= edges(b) & logx <= edges(b + 1));
        end
        if numel(idxb) < 3
          continue;
        end
        bx(end + 1, 1) = 10 ^ mean(logx(idxb)); %#ok<AGROW>
        by(end + 1, 1) = mean(my(idxb)); %#ok<AGROW>
      end
      if numel(bx) > 1
        [bx, ord] = sort(bx);
        by = by(ord);
        h_curve = semilogx(bx, by ./ 1e9, '-o', 'LineWidth', 2.0, ...
                           'MarkerSize', 5, 'Color', [0.85 0.10 0.10], ...
                           'MarkerFaceColor', 'none');
      end
    end
  end

  set(gca, 'XScale', 'log');
  xlabel('Plastic strain rate (1/s)');
  ylabel('Flow stress (GPa)');
  title(strrep(case_names{i}, '_', '\_'), 'Interpreter', 'tex');
  grid on;
  box on;
  if ~isempty(h_scatter) && ~isempty(h_curve)
    legend([h_scatter, h_curve], ...
           {'Interface particle samples', 'Mean flow stress by strain-rate bin'}, ...
           'Location', 'northwest');
  elseif ~isempty(h_scatter)
    legend(h_scatter, {'Interface particle samples'}, 'Location', 'northwest');
  elseif ~isempty(h_curve)
    legend(h_curve, {'Mean flow stress by strain-rate bin'}, 'Location', 'northwest');
  end
  print(fig, out_png, '-dpng', '-r220');
  close(fig);
end
