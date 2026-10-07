function model = prepare_baseband_model_5L(Op,Od,Oc,Olo,g,gamma_t,D_RF, ...
    N_eff,mu12,eps0,h_bar,lam_p,lam_c,lam_d,T_atom,Nv)

m_Cs = 2.2069e-25;
k_B = 1.380649e-23;
u = sqrt(2*k_B*T_atom/m_Cs);
v_array = linspace(-3*u,3*u,Nv);
dv = v_array(2)-v_array(1);
P_v = exp(-(v_array/u).^2)/(sqrt(pi)*u);
weights = P_v*dv;

kp = 2*pi/lam_p;
kd = 2*pi/lam_d;
kc = 2*pi/lam_c;
Gamma_pop = [0;g(1);g(2);g(3);g(4)];
Gamma_coh = build_coherence_decay(Gamma_pop,50e3*2*pi);

n = 5;
N = n^2;
trace_eq = zeros(1,N);
trace_eq([1,7,13,19,25]) = 1;
dL = rf_drive_liouvillian_derivative_5L();

L0 = cell(1,Nv);
rhs = cell(1,Nv);
chi_prime_dc = 0;
K = 2*N_eff*mu12^2/(eps0*h_bar*Op);

for iv = 1:Nv
    v = v_array(iv);
    D_vec = make_detuning_vector_5L(v,kp,kd,kc,D_RF);
    [~,rho0] = solve_steady_state_5L(Op,Od,Oc,Olo,Gamma_pop,Gamma_coh,gamma_t,D_vec);
    L0{iv} = build_liouvillian_5L(Op,Od,Oc,Olo,Gamma_pop,Gamma_coh,gamma_t,D_vec);
    rhs{iv} = dL*rho0;
    rhs{iv}(N) = 0;

    A0 = -L0{iv};
    A0(N,:) = trace_eq;
    drho0 = A0\rhs{iv};
    chi_prime_dc = chi_prime_dc + K*drho0(2)*weights(iv);
end

if ~isfinite(abs(chi_prime_dc)) || abs(chi_prime_dc) <= realmin
    error('RAQR:InvalidDCTransconductance', ...
        '5L DC transconductance is zero or non-finite; BCOD quadrature is undefined.');
end

psi = angle(chi_prime_dc);
[c_re,c_im] = probe_coherence_readout_vectors(n);
c_det = cos(psi)*c_re + sin(psi)*c_im;

model = struct();
model.levels = n;
model.L0 = L0;
model.rhs = rhs;
model.weights = weights;
model.trace_eq = trace_eq;
model.K = K;
model.c_det = c_det;
model.readout_phase_rad = psi;
model.chi_prime_dc = chi_prime_dc;
model.Nv = Nv;
end
