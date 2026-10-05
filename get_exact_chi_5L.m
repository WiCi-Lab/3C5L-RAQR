function [chi_0, chi_p, rho11_avg] = get_exact_chi_5L(Op, Od, Oc, Olo, g, gamma_t, D_RF, N_eff, mu12, eps0, h_bar, lam_p, lam_c, lam_d, T_BB, Nv)
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

d_Olo = 1e-4 * (abs(Olo) + 1e-6);
rho21_0_array = zeros(1, Nv);
rho21_p_array = zeros(1, Nv);
rho21_m_array = zeros(1, Nv);
rho11_array = zeros(1, Nv);

for idx = 1:Nv
    v = v_array(idx); weight = P_v(idx) * dv;
    D_vec = make_detuning_vector_5L(v,kp,kd,kc,D_RF);

    [~, rho_ss_0] = solve_steady_state_5L(Op, Od, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec);
    rho21_0_array(idx) = rho_ss_0(2) * weight;
    rho11_array(idx)   = real(rho_ss_0(1)) * weight;

    [~, rho_ss_p] = solve_steady_state_5L(Op, Od, Oc, Olo + d_Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec);
    [~, rho_ss_m] = solve_steady_state_5L(Op, Od, Oc, Olo - d_Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec);
    rho21_p_array(idx) = rho_ss_p(2) * weight;
    rho21_m_array(idx) = rho_ss_m(2) * weight;
end

K = (N_eff * mu12^2) / (eps0 * h_bar * Op);
chi_0 = K * sum(rho21_0_array);
chi_p = K * (sum(rho21_p_array) - sum(rho21_m_array)) / (2*d_Olo);
rho11_avg = sum(rho11_array);
end
