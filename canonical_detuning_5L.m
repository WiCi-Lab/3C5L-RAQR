function d = canonical_detuning_5L(detuning)
if isscalar(detuning)
    d=[0,0,0,detuning];
elseif numel(detuning)==4
    d=reshape(detuning,1,[]);
else
    error('5L detuning must be scalar or [Delta_p, Delta_d, Delta_c, Delta_RF].');
end
end
