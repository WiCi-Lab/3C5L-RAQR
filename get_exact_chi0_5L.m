function chi_0 = get_exact_chi0_5L(Op, Od, Oc, Olo, g, gamma_t, D_RF, ...
    N_eff, mu12, eps0, h_bar, lam_p, lam_c, lam_d, T_atom, Nv)
m_Cs = 2.2069e-25; k_B = 1.380649e-23;
u = sqrt(2*k_B*T_atom/m_Cs);
v_array = linspace(-3*u, 3*u, Nv);
dv = v_array(2)-v_array(1);
P_v = exp(-(v_array/u).^2)/(sqrt(pi)*u);

kp = 2*pi/lam_p; kd = 2*pi/lam_d; kc = 2*pi/lam_c;
Gamma_pop = [0; g(1); g(2); g(3); g(4)];
Gamma_coh = build_coherence_decay(Gamma_pop, 50e3*2*pi);
rho21_sum = 0;
for idx = 1:Nv
    v = v_array(idx);
    D_vec = make_detuning_vector_5L(v,kp,kd,kc,D_RF);
    [~, rho_ss] = solve_steady_state_5L(Op, Od, Oc, Olo, Gamma_pop, ...
        Gamma_coh, gamma_t, D_vec);
    rho21_sum = rho21_sum + rho_ss(2)*P_v(idx)*dv;
end
K = 2*N_eff*mu12^2/(eps0*h_bar*Op);
chi_0 = K*rho21_sum;
end
