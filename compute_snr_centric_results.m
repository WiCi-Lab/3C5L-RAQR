function r = compute_snr_centric_results(p,cfg)

cfg_snr = cfg;
cfg_snr.compute_communication_bandwidth = false;

P_t_W = 10.^((p.P_t_dBm-30)/10);
E_inc = sqrt(2*p.Z_0*P_t_W*p.G_t/(4*pi*p.dist^2));

[SNR_4L,~,~,~,~,~,~,EIT_FWHM_4L,SNR_4L_fixedLO] = ...
    get_system_metrics('4L',E_inc,p.h_bar,p.epsilon_0,p.Z_0,p.q,p.k_B,p.T_BB,p.L_cell,p.N_eff_4L,p.Upsilon_4L,p.P_p_4L,p.P_c_4L, ...
    p.r_beam_4L,p.E_LO_RF_4L,p.mu_12_4L,p.mu_23_4L,p.mu_RF_4L, ...
    p.lam_p_4L,p.lam_c_4L,p.gammas_4L_base,p.G_amp_V,p.alpha_pd, ...
    p.P_LO_opt,p.R_load,p.T_sys_PD,p.B_Signal_Common,p.detuning_RF_4L, ...
    p.N_v_4L,cfg_snr);

[SNR_5L,~,~,~,~,~,~,EIT_FWHM_5L,SNR_5L_fixedLO] = ...
    get_system_metrics('5L',E_inc,p.h_bar,p.epsilon_0,p.Z_0,p.q,p.k_B,p.T_BB,p.L_cell,p.N_eff_5L,p.Upsilon_5L,p.P_p_5L,p.P_c_5L, ...
    p.r_beam_5L,p.E_LO_RF_5L,p.mu_12_5L,p.mu_23_5L,p.mu_RF_5L, ...
    p.lam_p_5L,p.lam_c_5L,p.gammas_5L_base,p.G_amp_V,p.alpha_pd, ...
    p.P_LO_opt,p.R_load,p.T_sys_PD,p.B_Signal_Common,p.detuning_RF_5L, ...
    p.N_v_5L,p.P_d_5L,p.mu_34_5L,p.lam_d_5L,cfg_snr);

[SNR_CL,~,~] = get_system_metrics('CL',E_inc,p.h_bar,p.epsilon_0, ...
    p.Z_0,p.q,p.k_B,p.T_BB,p.L_cell,p.N_0_Base, ...
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
    p.B_Signal_Common,0,0,cfg_snr);

linear4 = SNR_4L_fixedLO(1).*10.^((p.P_t_dBm-p.P_t_dBm(1))/10);
linear5 = SNR_5L_fixedLO(1).*10.^((p.P_t_dBm-p.P_t_dBm(1))/10);
comp4 = 10*log10(max(SNR_4L_fixedLO./linear4,realmin));
comp5 = 10*log10(max(SNR_5L_fixedLO./linear5,realmin));
idx4 = find(comp4<=-cfg.compression_point_db,1,'first');
idx5 = find(comp5<=-cfg.compression_point_db,1,'first');
P1dB_4L = NaN; P1dB_5L = NaN;
if ~isempty(idx4), P1dB_4L = p.P_t_dBm(idx4); end
if ~isempty(idx5), P1dB_5L = p.P_t_dBm(idx5); end

r = struct();
r.P_t_dBm = p.P_t_dBm;
r.SNR_CL = SNR_CL;
r.SNR_4L = SNR_4L;
r.SNR_5L = SNR_5L;
r.SNR_4L_fixedLO = SNR_4L_fixedLO;
r.SNR_5L_fixedLO = SNR_5L_fixedLO;
r.EIT_FWHM_4L = EIT_FWHM_4L;
r.EIT_FWHM_5L = EIT_FWHM_5L;
r.P1dB_4L = P1dB_4L;
r.P1dB_5L = P1dB_5L;
end
