function cfg = default_model_cfg()
cfg = struct();

cfg.raqr_nonlinearity = 'liouvillian_rapp';
cfg.linear_fraction = 0.1;
cfg.nonlinear_phase_samples = 32;
cfg.nonlinear_velocity_points = 61;
cfg.compression_point_db = 1;
cfg.rapp_smoothness = 2.5;
cfg.compression_fallback_fraction = 0.5;

cfg.eit_linewidth_scan_half_span_Hz = 20e6;
cfg.eit_linewidth_scan_points = 121;
cfg.eit_linewidth_velocity_points = 61;
cfg.eit_linewidth_max_half_span_Hz = 320e6;
cfg.eit_linewidth_max_expansions = 5;
cfg.eit_linewidth_min_relative_contrast = 1e-6;

cfg.compute_communication_bandwidth = true;
cfg.communication_bandwidth_method = 'lo_biased_3db';
cfg.communication_bandwidth_frequency_grid_Hz = ...
    [0, logspace(2,log10(100e6),140)];
cfg.communication_bandwidth_velocity_points = 101;
cfg.communication_bandwidth_refine_tolerance = 1e-3;
cfg.communication_bandwidth_refine_max_iterations = 30;
end
