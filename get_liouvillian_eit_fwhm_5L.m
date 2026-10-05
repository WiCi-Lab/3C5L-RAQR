function BW_Hz = get_liouvillian_eit_fwhm_5L(Op,Od,Oc,Olo,g,gamma_t, ...
    detuning,N_eff,mu12,eps0,h_bar,lam_p,lam_c,lam_d,T_atom,Nv,cfg)

persistent cache5
if isempty(cache5)
    cache5=containers.Map('KeyType','char','ValueType','double');
end

d0=canonical_detuning_5L(detuning);
Nf=make_odd_grid_size(cfg.eit_linewidth_scan_points,41);
vals=[Op,Od,Oc,Olo,reshape(g,1,[]),gamma_t,d0,N_eff,mu12,eps0,h_bar, ...
    lam_p,lam_c,lam_d,T_atom,Nv,cfg.eit_linewidth_scan_half_span_Hz,Nf, ...
    cfg.eit_linewidth_max_half_span_Hz,cfg.eit_linewidth_max_expansions, ...
    cfg.eit_linewidth_min_relative_contrast];
key=make_eit_cache_key('5L',vals);
if isKey(cache5,key)
    BW_Hz=cache5(key);
    return;
end

half_span=min(max(cfg.eit_linewidth_scan_half_span_Hz,1e3), ...
    cfg.eit_linewidth_max_half_span_Hz);
BW_Hz=NaN;
for iexp=1:max(1,round(cfg.eit_linewidth_max_expansions))
    f_scan=linspace(-half_span,half_span,Nf);
    absorption=zeros(size(f_scan));
    for kf=1:numel(f_scan)
        dk=d0;
        dk(3)=dk(3)+2*pi*f_scan(kf);
        chi=get_exact_chi0_5L(Op,Od,Oc,Olo,g,gamma_t,dk,N_eff,mu12, ...
            eps0,h_bar,lam_p,lam_c,lam_d,T_atom,Nv);
        absorption(kf)=-imag(chi);
    end
    [BW_try,ok]=extract_main_feature_fwhm(f_scan,absorption, ...
        cfg.eit_linewidth_min_relative_contrast);
    if ok && BW_try<1.2*half_span
        BW_Hz=BW_try;
        break;
    end
    if half_span>=cfg.eit_linewidth_max_half_span_Hz
        break;
    end
    half_span=min(2*half_span,cfg.eit_linewidth_max_half_span_Hz);
end

if ~isfinite(BW_Hz) || BW_Hz<=0
    error('RAQR:UnresolvedEITLinewidth', ...
        ['3C5L full-Liouvillian spectrum did not yield a resolved EIT/EIA FWHM. ' ...
         'Increase the linewidth scan span/points or velocity grid.']);
end
cache5(key)=BW_Hz;
end
