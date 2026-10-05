function [BW_Hz,H_grid,info] = find_3db_bandwidth_adaptive(model,f_grid,tol_rel,max_iter)

if nargin < 3 || isempty(tol_rel), tol_rel = 1e-3; end
if nargin < 4 || isempty(max_iter), max_iter = 30; end

f = abs(f_grid(:));
f = unique(sort(f),'stable');
if isempty(f) || f(1) ~= 0
    f = [0;f];
end

H_grid = evaluate_baseband_model(model,f);
mag = abs(H_grid(:));
if ~isfinite(mag(1)) || mag(1) <= realmin
    BW_Hz = NaN;
    info = struct('iterations',0,'bracket',[NaN NaN],'target',NaN);
    warning('RAQR:InvalidBasebandResponse', ...
        'Cannot determine a DC-normalized -3 dB bandwidth.');
    return;
end

target = mag(1)/sqrt(2);
g = mag-target;
idx = find(g <= 0,1,'first');
if isempty(idx)
    BW_Hz = NaN;
    info = struct('iterations',0,'bracket',[f(end) NaN],'target',target);
    warning('RAQR:BandwidthOutsideGrid', ...
        'The -3 dB crossing lies above the simulated frequency grid.');
    return;
elseif idx == 1
    BW_Hz = 0;
    info = struct('iterations',0,'bracket',[0 0],'target',target);
    return;
end

flo = f(idx-1);
fhi = f(idx);
glo = g(idx-1);
ghi = g(idx);
initial_bracket = [flo fhi];

for it = 1:max_iter
    if flo > 0
        fmid = sqrt(flo*fhi);
    else
        fmid = (flo+fhi)/2;
    end
    Hmid = evaluate_baseband_model(model,fmid);
    gmid = abs(Hmid)-target;

    if gmid > 0
        flo = fmid;
        glo = gmid;
    else
        fhi = fmid;
        ghi = gmid;
    end

    if (fhi-flo)/max(fhi,realmin) <= tol_rel
        break;
    end
end

if flo > 0 && isfinite(glo) && isfinite(ghi) && ghi ~= glo
    xlo = log10(flo); xhi = log10(fhi);
    x = xlo + (0-glo)*(xhi-xlo)/(ghi-glo);
    BW_Hz = 10^x;
else
    BW_Hz = (flo+fhi)/2;
end

info = struct();
info.iterations = it;
info.initial_bracket = initial_bracket;
info.final_bracket = [flo fhi];
info.target = target;
info.relative_bracket_width = (fhi-flo)/max(fhi,realmin);
end
