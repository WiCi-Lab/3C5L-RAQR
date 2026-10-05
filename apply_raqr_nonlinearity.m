function [E_eff, power_ratio, fixed_lo_ratio] = apply_raqr_nonlinearity(E_inc, E_LO_RF, muRF, ...
    h_bar, linewidth_Hz, P_probe, P_out_LO, chi_LO, chi_prime_LO, ...
    lam_p, L_cell, chi_solver, model_cfg)

switch lower(model_cfg.raqr_nonlinearity)
    case 'liouvillian_rapp'
        E_linewidth = h_bar * (2*pi*linewidth_Hz) / muRF;
        E_eff = E_inc;
        fixed_lo_ratio = exact_harmonic_power_ratio(E_inc, E_LO_RF, muRF, h_bar, ...
            E_linewidth, P_probe, P_out_LO, chi_LO, chi_prime_LO, lam_p, ...
            L_cell, chi_solver, model_cfg);
        power_ratio = liouvillian_calibrated_rapp(E_inc, fixed_lo_ratio, ...
            E_LO_RF, E_linewidth, model_cfg);

    case 'linear'
        E_eff = E_inc;
        power_ratio = ones(size(E_inc));
        fixed_lo_ratio = power_ratio;

    otherwise
        error(['Unsupported RAQR nonlinearity in the clean release: %s. ' ...
            'Use liouvillian_rapp or linear.'], model_cfg.raqr_nonlinearity);
end
end
