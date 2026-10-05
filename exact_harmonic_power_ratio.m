function ratio = exact_harmonic_power_ratio(E_inc, E_LO_RF, muRF, h_bar, ...
    E_linewidth, P_probe, P_out_LO, chi_LO, chi_prime_LO, lam_p, ...
    L_cell, chi_solver, model_cfg)

ratio = ones(size(E_inc));
N_phase = max(8, round(model_cfg.nonlinear_phase_samples));
if mod(N_phase, 2) ~= 0
    N_phase = N_phase + 1;
end
theta = 2*pi*(0:N_phase-1)/N_phase;
c_theta = cos(theta);

phi_LO = pi * L_cell * real(chi_LO) / lam_p;
g_LO = sqrt(max(P_out_LO, 0)) * exp(1i*phi_LO);
dg_dOmega = g_LO * (pi*L_cell/lam_p) * ...
    (imag(chi_prime_LO) + 1i*real(chi_prime_LO));
dg_mag = abs(dg_dOmega);
if dg_mag <= realmin
    ratio(:) = 0;
    return;
end
detector_phase = angle(dg_dOmega);

E_small = model_cfg.linear_fraction * min(abs(E_LO_RF), E_linewidth);

for idx = 1:numel(E_inc)
    E_sig = abs(E_inc(idx));
    if E_sig == 0 || E_sig <= E_small
        ratio(idx) = 1;
        continue;
    end

    Omega_sig = muRF * E_sig / h_bar;
    q_theta = zeros(1, N_phase);
    for k = 1:N_phase
        E_total = E_LO_RF + E_sig*c_theta(k);
        Omega_total = muRF * E_total / h_bar;
        chi_k = chi_solver(Omega_total);

        log_amp = pi * L_cell * min(imag(chi_k), 0) / lam_p;
        optical_field = sqrt(P_probe) * exp(log_amp) * ...
            exp(1i*pi*L_cell*real(chi_k)/lam_p);
        q_theta(k) = real(optical_field * exp(-1i*detector_phase));
    end

    harmonic_amp = (2/N_phase) * sum(q_theta .* c_theta);
    linear_amp = Omega_sig * dg_mag;
    ratio(idx) = max(real((harmonic_amp/linear_amp)^2), 0);
end
end
