function [BW_Hz, tau_f] = get_exact_transient_bandwidth_5L(Op, Od, Oc, Olo, g, gamma_t, D_RF, T_BB, lam_p, lam_c, lam_d, Nv)
m_Cs = 2.2069e-25; k_B = 1.380649e-23;
u = sqrt(2 * k_B * T_BB / m_Cs);
v_array = linspace(-3*u, 3*u, Nv); dv = v_array(2) - v_array(1);
P_v = (1 / (sqrt(pi) * u)) * exp(-(v_array / u).^2);

kp = 2*pi/lam_p; kd = 2*pi/lam_d; kc = 2*pi/lam_c;
Gamma_pop = [0; g(1); g(2); g(3); g(4)];
Gamma_coh = zeros(5, 5);
for i = 1:5
    for j = 1:5
        Gamma_coh(i,j) = (Gamma_pop(i) + Gamma_pop(j))/2;
        if i ~= j
            Gamma_coh(i,j) = Gamma_coh(i,j) + 50e3*2*pi;
        end
    end
end

tau_physical_limit = 1 / gamma_t;
t_max = 3 * tau_physical_limit;

t_out = linspace(0, t_max, 500)';
im_rho12_macro = zeros(length(t_out), 1); eye25 = eye(25);

for idx = 1:Nv
    v = v_array(idx); weight = P_v(idx) * dv;
    D_vec = make_detuning_vector_5L(v,kp,kd,kc,D_RF);
    [~, rho_ss_ON] = solve_steady_state_5L(Op, Od, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec);

    L_off = zeros(25, 25);
    for k = 1:25, L_off(:, k) = bloch_5L_matrix(0, eye25(:, k), Op, Od, Oc, 0, Gamma_pop, Gamma_coh, gamma_t, D_vec); end

    [V, D_eig] = eig(L_off);
    invV_rho0 = V \ rho_ss_ON;
    exp_Dt = exp(diag(D_eig) * t_out');
    rho_t = V * (exp_Dt .* invV_rho0);
    im_rho12_macro = im_rho12_macro + (-imag(rho_t(2, :).')) * weight;
end

val_start = im_rho12_macro(1); val_end = im_rho12_macro(end);
val_target = val_end + (val_start - val_end) * exp(-1);

if abs(val_start - val_end) < 1e-15
    tau_f = 1 / (max(Gamma_pop(2:end)) + gamma_t);
else
    idx_cross = find(ifelse(val_start > val_end, im_rho12_macro <= val_target, im_rho12_macro >= val_target), 1, 'first');
    if isempty(idx_cross) || idx_cross == 1
        tau_f = t_out(end);
    else
        t1 = t_out(idx_cross-1); t2 = t_out(idx_cross);
        v1 = im_rho12_macro(idx_cross-1); v2 = im_rho12_macro(idx_cross);
        tau_f = t1 + (val_target - v1) * (t2 - t1) / (v2 - v1);
    end
end
BW_Hz = 1 / (2 * pi * tau_f);
end
