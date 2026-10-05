clc;
clear;
close all;
rng(42,'twister');

p = parameters();
cfg = config();

r = compute_snr_centric_results(p,cfg);

P_t_dBm = r.P_t_dBm;
SNR_CL = r.SNR_CL;
SNR_4L = r.SNR_4L;
SNR_5L = r.SNR_5L;
SNR_4L_fixedLO = r.SNR_4L_fixedLO;
SNR_5L_fixedLO = r.SNR_5L_fixedLO;
P1dB_4L = r.P1dB_4L;
P1dB_5L = r.P1dB_5L;

figure;
plot(P_t_dBm,10*log10(max(SNR_CL,realmin)),'k-o','LineWidth',2, ...
    'DisplayName','Classical RF');
hold on;
plot(P_t_dBm,10*log10(max(SNR_4L,realmin)),'b-s','LineWidth',2, ...
    'DisplayName','2C4L-RAQR');
plot(P_t_dBm,10*log10(max(SNR_4L_fixedLO,realmin)),'b-.s', ...
    'LineWidth',2,'MarkerFaceColor','w','HandleVisibility','off');
plot(P_t_dBm,10*log10(max(SNR_5L,realmin)),'r-d','LineWidth',2, ...
    'DisplayName','3C5L-RAQR');
plot(P_t_dBm,10*log10(max(SNR_5L_fixedLO,realmin)),'r--d', ...
    'LineWidth',2,'MarkerFaceColor','w','HandleVisibility','off');

if isfinite(P1dB_4L)
    xline(P1dB_4L,'b:','2C4L 1-dB','LineWidth',1.2, ...
        'LabelVerticalAlignment','bottom','HandleVisibility','off');
end
if isfinite(P1dB_5L)
    xline(P1dB_5L,'r:','3C5L 1-dB','LineWidth',1.2, ...
        'LabelVerticalAlignment','top','HandleVisibility','off');
end

xlabel('Transmit Power $P_t$ (dBm)','Interpreter','latex','FontSize',12);
ylabel('SNR (dB)','Interpreter','latex','FontSize',12);
legend('Location','best','FontSize',12,'Interpreter','latex');
grid on;
box on;
set(gca,'FontSize',11,'XMinorGrid','on','YMinorGrid','on');

ax_ins = axes('Position',[0.17,0.65,0.24,0.24]);
hold(ax_ins,'on');
plot(ax_ins,P_t_dBm,10*log10(max(SNR_CL,realmin)),'k-o', ...
    'LineWidth',1.5,'MarkerSize',4);
plot(ax_ins,P_t_dBm,10*log10(max(SNR_4L,realmin)),'b-s', ...
    'LineWidth',1.5,'MarkerSize',4);
plot(ax_ins,P_t_dBm,10*log10(max(SNR_4L_fixedLO,realmin)),'b--s', ...
    'LineWidth',1.2,'MarkerSize',4,'MarkerFaceColor','w');
xlim(ax_ins,[-40,-35]);
grid(ax_ins,'on');
box(ax_ins,'on');
set(ax_ins,'FontSize',8);
