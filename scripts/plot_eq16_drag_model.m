function plot_eq16_drag_model()
  script_dir = fileparts(mfilename("fullpath"));
  root_dir = fileparts(script_dir);
  input_json = fullfile(root_dir, "examples", "mechanics", "inputs", "simple_impact_thermal.json");
  fig_dir = fullfile(root_dir, "latex", "figures");

  if exist(fig_dir, "dir") ~= 7
    mkdir(fig_dir);
  end

  inputs = jsondecode(fileread(input_json));
  materials = {
    build_material(inputs, 1, "Al projectile"),
    build_material(inputs, 2, "Cu substrate")
  };

  epsdot_vals = logspace(3, 8, 250);
  eps_p_fixed = 0.20;
  eps_p_vals = linspace(0.0, 0.50, 250);
  strain_curve_rates = [1.0e4, 1.0e6, 1.0e7];

  rate_fig = figure("visible", "off", "position", [100, 100, 1200, 480]);
  for i = 1:numel(materials)
    mat = materials{i};
    subplot(1, 2, i);
    sigma_low = arrayfun(@(r) jc_low_stress(mat, eps_p_fixed, r, mat.T_ref), epsdot_vals);
    sigma_total = arrayfun(@(r) jc_eq16_stress(mat, eps_p_fixed, r, mat.T_ref), epsdot_vals);
    loglog(epsdot_vals, sigma_low * 1.0e-6, "--", "linewidth", 2.0);
    hold on;
    loglog(epsdot_vals, sigma_total * 1.0e-6, "-", "linewidth", 2.0);
    grid on;
    xlabel("Plastic strain rate, \epsilon_p dot (s^{-1})");
    ylabel("Flow stress (MPa)");
    title(sprintf("%s, \\epsilon_p = %.2f, T = %.0f K", mat.name, eps_p_fixed, mat.T_ref));
    legend("Low-rate JC", "JC + Eq.16 drag", "location", "northwest");
  end
  set(rate_fig, "paperpositionmode", "auto");
  print(rate_fig, fullfile(fig_dir, "jc_eq16_stress_vs_strain_rate.png"), "-dpng", "-r200");
  print(rate_fig, fullfile(fig_dir, "jc_eq16_stress_vs_strain_rate.pdf"), "-dpdf", "-bestfit");
  close(rate_fig);

  strain_fig = figure("visible", "off", "position", [100, 100, 1200, 480]);
  for i = 1:numel(materials)
    mat = materials{i};
    subplot(1, 2, i);
    hold on;
    for rate = strain_curve_rates
      sigma_total = arrayfun(@(eps_p) jc_eq16_stress(mat, eps_p, rate, mat.T_ref), eps_p_vals);
      plot(eps_p_vals, sigma_total * 1.0e-6, "linewidth", 2.0, ...
           "displayname", sprintf("\\epsilon_p dot = %.1e s^{-1}", rate));
    end
    grid on;
    xlabel("Equivalent plastic strain, \epsilon_p");
    ylabel("Flow stress (MPa)");
    title(sprintf("%s, T = %.0f K", mat.name, mat.T_ref));
    legend("location", "northwest");
  end
  set(strain_fig, "paperpositionmode", "auto");
  print(strain_fig, fullfile(fig_dir, "jc_eq16_stress_vs_strain.png"), "-dpng", "-r200");
  print(strain_fig, fullfile(fig_dir, "jc_eq16_stress_vs_strain.pdf"), "-dpdf", "-bestfit");
  close(strain_fig);

  verification_path = fullfile(fig_dir, "jc_eq16_drag_verification.tsv");
  fid = fopen(verification_path, "w");
  fprintf(fid, "material\tdrag_slope_Pa_s\tBd_Pa_s\trho_mobile_m^-2\tburgers_m\tsigma_drag_1e7_MPa\tratio_drag_2x\tclosure_residual_Pa\n");
  for i = 1:numel(materials)
    mat = materials{i};
    sample_rate = 1.0e7;
    sigma_low = jc_low_stress(mat, eps_p_fixed, sample_rate, mat.T_ref);
    sigma_drag_1 = drag_stress(mat, sample_rate);
    sigma_drag_2 = drag_stress(mat, 2.0 * sample_rate);
    sigma_total = jc_eq16_stress(mat, eps_p_fixed, sample_rate, mat.T_ref);
    slope = mat.drag_Bd / (mat.drag_rho_mobile * mat.drag_burgers ^ 2);
    residual = sigma_total - sigma_low - sigma_drag_1;
    fprintf(fid, "%s\t%.8e\t%.8e\t%.8e\t%.8e\t%.8f\t%.8f\t%.8e\n", ...
            mat.name, slope, mat.drag_Bd, mat.drag_rho_mobile, mat.drag_burgers, ...
            sigma_drag_1 * 1.0e-6, sigma_drag_2 / sigma_drag_1, residual);
  end
  fclose(fid);
end

function mat = build_material(inputs, idx, name)
  mat.name = name;
  mat.A = inputs.yield_stress.value(idx);
  mat.B = inputs.jc_B.value(idx);
  mat.n = inputs.jc_n.value(idx);
  mat.C = inputs.jc_C.value(idx);
  mat.epsdot0 = inputs.jc_epsdot0.value(idx);
  mat.T_ref = inputs.reference_temperature.value(idx);
  mat.T_melt = inputs.jc_Tmelt.value(idx);
  mat.m_thermal = inputs.jc_m.value(idx);
  mat.drag_Bd = inputs.drag_Bd.value(idx);
  mat.drag_rho_mobile = inputs.drag_mobile_dislocation_density.value(idx);
  mat.drag_burgers = inputs.drag_burgers_vector.value(idx);
end

function sigma = jc_eq16_stress(mat, eps_p, eps_p_dot, T)
  sigma = jc_low_stress(mat, eps_p, eps_p_dot, T) + drag_stress(mat, eps_p_dot);
end

function sigma = jc_low_stress(mat, eps_p, eps_p_dot, T)
  hardening = mat.A + mat.B * (eps_p ^ mat.n);
  if mat.C == 0.0 || mat.epsdot0 <= 0.0
    rate_factor = 1.0;
  else
    eps_eff = max(eps_p_dot, mat.epsdot0);
    rate_factor = 1.0 + mat.C * log(eps_eff / mat.epsdot0);
  end

  if mat.m_thermal <= 0.0 || mat.T_melt <= mat.T_ref
    thermal_factor = 1.0;
  else
    T_star = (T - mat.T_ref) / (mat.T_melt - mat.T_ref);
    T_star = min(max(T_star, 0.0), 1.0);
    thermal_factor = max(0.0, 1.0 - T_star ^ mat.m_thermal);
  end

  sigma = hardening * rate_factor * thermal_factor;
end

function sigma = drag_stress(mat, eps_p_dot)
  if mat.drag_Bd <= 0.0 || mat.drag_rho_mobile <= 0.0 || mat.drag_burgers <= 0.0
    sigma = 0.0;
    return;
  end

  sigma = mat.drag_Bd * max(eps_p_dot, 0.0) / ...
          (mat.drag_rho_mobile * mat.drag_burgers ^ 2);
end

plot_eq16_drag_model();
