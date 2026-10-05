function ratio = liouvillian_calibrated_rapp(E_inc, raw_ratio, E_LO_RF, E_linewidth, cfg)
E = abs(E_inc(:)); raw = real(raw_ratio(:));
[E_sorted,order] = sort(E);
raw_sorted = raw(order);
if numel(raw_sorted) >= 3
    raw_sorted = movmedian(raw_sorted,3);
end
target = 10^(-cfg.compression_point_db/10);
valid_start = cfg.linear_fraction*min(abs(E_LO_RF),E_linewidth);
idx = find(E_sorted >= valid_start & raw_sorted <= target,1,'first');
if isempty(idx)
    E_comp = cfg.compression_fallback_fraction*min(abs(E_LO_RF),E_linewidth);
else
    E_comp = E_sorted(idx);
end
E_comp = max(E_comp,realmin);

p = cfg.rapp_smoothness;
E_sat = E_comp/(10^(cfg.compression_point_db*p/10)-1)^(1/(2*p));
ratio = (1+(abs(E_inc)/E_sat).^(2*p)).^(-1/p);
end
