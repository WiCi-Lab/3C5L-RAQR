function [BW_Hz,ok] = extract_main_feature_fwhm(f_Hz,spectrum,min_rel_contrast)
f=f_Hz(:); y=real(spectrum(:)); N=numel(f);
BW_Hz=NaN; ok=false;
if N<9 || any(~isfinite(f)) || any(~isfinite(y)) || any(diff(f)<=0)
    return;
end
nedge=max(3,min(floor((N-3)/4),round(0.08*N)));
bg_left=median(y(1:nedge));
bg_right=median(y(end-nedge+1:end));
background=bg_left+(bg_right-bg_left)*(f-f(1))/(f(end)-f(1));
amp=abs(y-background);
search_idx=(nedge+1):(N-nedge);
[peak_amp,iloc]=max(amp(search_idx));
ipeak=search_idx(iloc);
scale=max([max(abs(y)),abs(bg_left),abs(bg_right),realmin]);
if peak_amp<=max(min_rel_contrast,0)*scale
    return;
end
half_amp=peak_amp/2;
ileft=find(amp(1:ipeak)<=half_amp,1,'last');
iright_rel=find(amp(ipeak:end)<=half_amp,1,'first');
if isempty(ileft) || isempty(iright_rel)
    return;
end
iright=ipeak+iright_rel-1;
if ileft>=ipeak || iright<=ipeak
    return;
end
fleft=linear_level_crossing(f(ileft),f(ileft+1), ...
    amp(ileft),amp(ileft+1),half_amp);
fright=linear_level_crossing(f(iright-1),f(iright), ...
    amp(iright-1),amp(iright),half_amp);
BW_Hz=fright-fleft;
ok=isfinite(BW_Hz) && BW_Hz>0;
end
