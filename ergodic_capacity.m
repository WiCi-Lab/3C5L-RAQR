function cap = ergodic_capacity(SNR_vals, BW)

if isscalar(BW)
    BW_array = repmat(BW, size(SNR_vals));
else
    BW_array = BW;
end
cap = zeros(size(SNR_vals));
for i = 1:numel(SNR_vals)
    s = real(SNR_vals(i));
    if isnan(s) || isnan(BW_array(i))
        cap(i) = NaN;
    elseif s <= 0
        cap(i) = 0;
    elseif 1/s > 700
        cap(i) = BW_array(i) * s / log(2);
    else
        eta = exp(1/s) * expint(1/s) / log(2);
        cap(i) = BW_array(i) * real(eta);
    end
end
end
