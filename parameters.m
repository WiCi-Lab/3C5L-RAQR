function p = parameters()

p.h_bar = 1.0545718e-34;
p.epsilon_0 = 8.8541878e-12;
p.Z_0 = 376.73;
p.q = 1.60217663e-19;
p.a_0 = 5.2917721e-11;
p.k_B = 1.380649e-23;

p.T_BB = 290;
p.L_cell = 0.075;
p.N_0_Base = 1e10 * 1e7;
p.dist = 10;
p.G_t_dB = 10;
p.G_t = 10^(p.G_t_dB/10);
p.G_amp_dB = 30;
p.G_amp_V = 10^(p.G_amp_dB/20);
p.R_load = 1.0;
p.alpha_pd = 0.8;
p.P_LO_opt = 30e-3;
p.T_sys_PD = 290;
p.B_Signal_Common = 100e3;

p.Upsilon_4L = 0.02;
p.N_eff_4L = 1.5e10 * 1e7;
p.P_p_4L = 1e-3;
p.P_c_4L = 150e-3;
p.r_beam_4L = 0.38e-3;
p.E_LO_RF_4L = 0.03;
p.mu_12_4L = 2.59*p.q*p.a_0;
p.mu_23_4L = 0.022*p.q*p.a_0;
p.mu_RF_4L = 1443.48*p.q*p.a_0;
p.lam_p_4L = 852e-9;
p.lam_c_4L = 510e-9;
p.gammas_4L_base = [5.22e6*2*pi, 2.73e3*2*pi, 0.57e3*2*pi];
p.detuning_RF_4L = 0;

p.Upsilon_5L = 0.25;
p.N_eff_5L = 1.5e10 * 1e7;
p.P_p_5L = 50e-6;
p.P_d_5L = 20e-3;
p.P_c_5L = 40e-3;
p.r_beam_5L = 0.38e-3;
p.E_LO_RF_5L = 0.03;
p.mu_12_5L = 1.84*p.q*p.a_0;
p.mu_23_5L = 0.23*p.q*p.a_0;
p.mu_34_5L = 0.019*p.q*p.a_0;
p.mu_RF_5L = 1443.48*p.q*p.a_0;
p.lam_p_5L = 895e-9;
p.lam_d_5L = 636e-9;
p.lam_c_5L = 2245e-9;
p.gammas_5L_base = [4.56e6*2*pi, 0.98e6*2*pi, 0.57e3*2*pi, 2.73e3*2*pi];

p.detuning_RF_5L = -0.6e6*2*pi;

p.N_v_4L = 201;
p.N_v_5L = 201;
p.P_t_dBm = linspace(-80,20,25);
end
