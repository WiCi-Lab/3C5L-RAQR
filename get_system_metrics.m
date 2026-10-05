function [SNR, Cap, BW, P_out, SNR_PSN, SNR_QPN, rho11_avg, ...
    eit_linewidth, SNR_fixedLO, Cap_fixedLO] = get_system_metrics(type, ...
    E_inc_array, h_bar, epsilon_0, Z_0, q, k_B, T_BB, L_cell, N_eff_Base, ...
    Upsilon, P_p, P_c, r_beam, E_LO_RF, mu12, mu23, muRF, lam_p, lam_c, ...
    gammas_base, G_amp_V, alpha_pd, P_LO_opt, R_load, T_sys_PD, ...
    B_Signal_Common, opt_detuning_guess, N_v, varargin)

[model_cfg, varargin] = extract_model_cfg(varargin);
gamma_BBR_func = @(n_eff, T) bbr_decay_rate(n_eff, T);

if strcmp(type, '4L')
    gamma_BBR_47 = gamma_BBR_func(47.0, T_BB);
    gamma_BBR_48 = gamma_BBR_func(48.0, T_BB);
    gs = [gammas_base(1), gammas_base(2) + gamma_BBR_47, ...
        gammas_base(3) + gamma_BBR_48];

    Op = calculate_rabi(P_p, r_beam, mu12, Z_0, h_bar);
    Oc = calculate_rabi(P_c, r_beam, mu23, Z_0, h_bar);
    Olo = muRF * E_LO_RF / h_bar;

    m_Cs = 2.21e-25;
    v_th = sqrt(2 * k_B * T_BB / m_Cs);
    gamma_transit = v_th / r_beam;

    [chi_0, chi_p, rho11_avg] = get_exact_chi_4L(Op, Oc, Olo, gs, ...
        gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, T_BB, N_v);
    P_out = get_P_out_calculated(P_p, chi_0, lam_p, L_cell);

    Nv_eit = min(N_v, model_cfg.eit_linewidth_velocity_points);
    eit_linewidth = get_liouvillian_eit_fwhm_4L(Op, Oc, Olo, gs, ...
        gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, T_BB, Nv_eit, model_cfg);

    BW = select_communication_bandwidth_4L(model_cfg, N_v, Op, Oc, Olo, ...
        gs, gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, T_BB);
    if isscalar(BW), BW = repmat(BW, size(E_inc_array)); end

    chi_solver = @(Omega_RF) get_exact_chi0_4L(Op, Oc, Omega_RF, gs, ...
        gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, T_BB, model_cfg.nonlinear_velocity_points);
    [E_eff, nonlinear_ratio, fixed_lo_ratio] = apply_raqr_nonlinearity( ...
        E_inc_array, E_LO_RF, muRF, h_bar, eit_linewidth, P_p, P_out, ...
        chi_0, chi_p, lam_p, L_cell, chi_solver, model_cfg);

    mu_RF_val = muRF;
    r_b_val = r_beam;
    lam_p_val = lam_p;

elseif strcmp(type, '5L')
    if numel(varargin) < 3
        error('5L model requires P_d, mu34, and lam_d.');
    end
    P_d = varargin{1};
    mu34 = varargin{2};
    lam_d = varargin{3};

    gamma_BBR1 = gamma_BBR_func(47.0, T_BB);
    gamma_BBR2 = gamma_BBR_func(48.0, T_BB);
    gs = [gammas_base(1), gammas_base(2), ...
        gammas_base(3) + gamma_BBR2, gammas_base(4) + gamma_BBR1];

    Op = calculate_rabi(P_p, r_beam, mu12, Z_0, h_bar);
    Od = calculate_rabi(P_d, r_beam, mu23, Z_0, h_bar);
    Oc = calculate_rabi(P_c, r_beam, mu34, Z_0, h_bar);
    Olo = muRF * E_LO_RF / h_bar;

    m_Cs = 2.21e-25;
    v_th = sqrt(2 * k_B * T_BB / m_Cs);
    gamma_transit = v_th / r_beam;

    [chi_0, chi_p, rho11_avg] = get_exact_chi_5L(Op, Od, Oc, Olo, gs, ...
        gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, lam_d, T_BB, N_v);
    P_out = get_P_out_calculated(P_p, chi_0, lam_p, L_cell);

    Nv_eit = min(N_v, model_cfg.eit_linewidth_velocity_points);
    eit_linewidth = get_liouvillian_eit_fwhm_5L(Op, Od, Oc, Olo, gs, ...
        gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, lam_d, T_BB, Nv_eit, model_cfg);

    BW = select_communication_bandwidth_5L(model_cfg, N_v, Op, Od, Oc, Olo, ...
        gs, gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, lam_d, T_BB);
    if isscalar(BW), BW = repmat(BW, size(E_inc_array)); end

    chi_solver = @(Omega_RF) get_exact_chi0_5L(Op, Od, Oc, Omega_RF, gs, ...
        gamma_transit, opt_detuning_guess, N_eff_Base, mu12, epsilon_0, ...
        h_bar, lam_p, lam_c, lam_d, T_BB, model_cfg.nonlinear_velocity_points);
    [E_eff, nonlinear_ratio, fixed_lo_ratio] = apply_raqr_nonlinearity( ...
        E_inc_array, E_LO_RF, muRF, h_bar, eit_linewidth, P_p, P_out, ...
        chi_0, chi_p, lam_p, L_cell, chi_solver, model_cfg);

    mu_RF_val = muRF;
    r_b_val = r_beam;
    lam_p_val = lam_p;

elseif strcmp(type, 'CL')
    BW = 10e6 * ones(size(E_inc_array));
    A_eff = ((3e8/6.94e9)^2 / (4*pi)) * 10^(5.5/10);
    rho = A_eff * 10^(60/10);
    N_tot_val = k_B * T_BB * 10^(60/10) * 10^(6/10);
    P_sig = rho * E_inc_array.^2 / (2 * Z_0);
    SNR = P_sig / (N_tot_val * B_Signal_Common);
    Cap = ergodic_capacity(P_sig ./ (N_tot_val .* BW), BW);
    P_out = P_sig;
    SNR_PSN = SNR;
    SNR_QPN = SNR;
    SNR_fixedLO = SNR;
    Cap_fixedLO = Cap;
    rho11_avg = NaN;
    eit_linewidth = NaN;
    return;
else
    error('Unknown receiver type: %s', type);
end

kappa_2 = (pi * L_cell * mu_RF_val / (lam_p_val * h_bar)) * abs(chi_p);
rho = 4 * (G_amp_V * alpha_pd)^2 * Z_0 * P_LO_opt * P_out * kappa_2^2;

if Upsilon <= 0 || Upsilon > 1
    error('Upsilon must satisfy 0 < Upsilon <= 1.');
end
N_atoms = Upsilon * N_eff_Base * (pi * r_b_val^2 * L_cell);

N_thermal = k_B * T_sys_PD * G_amp_V^2;
N_shot = 2*q*R_load*G_amp_V^2*alpha_pd*(P_LO_opt + P_out);
N_qpn = (rho/(2*Z_0)) * (h_bar/mu_RF_val)^2 * ...
    (eit_linewidth*2*pi/N_atoms);
N_tot = N_thermal + N_shot + N_qpn;

P_sig = (rho * E_eff.^2 / (2 * Z_0)) .* nonlinear_ratio;
P_sig_fixedLO = (rho * E_inc_array.^2 / (2 * Z_0)) .* fixed_lo_ratio;

SNR_PSN = P_sig ./ (N_shot .* B_Signal_Common);
SNR_QPN = P_sig ./ (N_qpn .* B_Signal_Common);
SNR = P_sig ./ (N_tot .* B_Signal_Common);
Cap = ergodic_capacity(P_sig ./ (N_tot .* BW), BW);
SNR_fixedLO = P_sig_fixedLO ./ (N_tot .* B_Signal_Common);
Cap_fixedLO = ergodic_capacity(P_sig_fixedLO ./ (N_tot .* BW), BW);
end

function BW = select_communication_bandwidth_4L(cfg, N_v, Op, Oc, Olo, gs, ...
    gamma_transit, detuning, N_eff, mu12, eps0, h_bar, lam_p, lam_c, T_atom)
if ~cfg.compute_communication_bandwidth
    BW = NaN;
    return;
end
switch lower(cfg.communication_bandwidth_method)
    case 'lo_biased_3db'
        f_bw = cfg.communication_bandwidth_frequency_grid_Hz;
        Nv_bw = min(N_v, cfg.communication_bandwidth_velocity_points);
        model_bw = prepare_baseband_model_4L(Op, Oc, Olo, gs, gamma_transit, ...
            detuning, N_eff, mu12, eps0, h_bar, lam_p, lam_c, T_atom, Nv_bw);
        BW = find_3db_bandwidth_adaptive(model_bw, f_bw, ...
            cfg.communication_bandwidth_refine_tolerance, ...
            cfg.communication_bandwidth_refine_max_iterations);
        assert_valid_communication_bandwidth(BW, f_bw, '2C4L');
    case 'rf_off_1e'
        [BW,~] = get_exact_transient_bandwidth_4L(Op, Oc, Olo, gs, ...
            gamma_transit, detuning, T_atom, lam_p, lam_c, N_v);
    otherwise
        error(['Unsupported bandwidth method in the clean release: %s. ' ...
            'Use lo_biased_3db or rf_off_1e.'], cfg.communication_bandwidth_method);
end
end

function BW = select_communication_bandwidth_5L(cfg, N_v, Op, Od, Oc, Olo, gs, ...
    gamma_transit, detuning, N_eff, mu12, eps0, h_bar, lam_p, lam_c, lam_d, T_atom)
if ~cfg.compute_communication_bandwidth
    BW = NaN;
    return;
end
switch lower(cfg.communication_bandwidth_method)
    case 'lo_biased_3db'
        f_bw = cfg.communication_bandwidth_frequency_grid_Hz;
        Nv_bw = min(N_v, cfg.communication_bandwidth_velocity_points);
        model_bw = prepare_baseband_model_5L(Op, Od, Oc, Olo, gs, gamma_transit, ...
            detuning, N_eff, mu12, eps0, h_bar, lam_p, lam_c, lam_d, T_atom, Nv_bw);
        BW = find_3db_bandwidth_adaptive(model_bw, f_bw, ...
            cfg.communication_bandwidth_refine_tolerance, ...
            cfg.communication_bandwidth_refine_max_iterations);
        assert_valid_communication_bandwidth(BW, f_bw, '3C5L');
    case 'rf_off_1e'
        [BW,~] = get_exact_transient_bandwidth_5L(Op, Od, Oc, Olo, gs, ...
            gamma_transit, detuning, T_atom, lam_p, lam_c, lam_d, N_v);
    otherwise
        error(['Unsupported bandwidth method in the clean release: %s. ' ...
            'Use lo_biased_3db or rf_off_1e.'], cfg.communication_bandwidth_method);
end
end
