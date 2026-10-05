function d = canonical_detuning_4L(detuning)
if isscalar(detuning)
    d=[0,0,detuning];
elseif numel(detuning)==3
    d=reshape(detuning,1,[]);
else
    error('4L detuning must be scalar or [Delta_p, Delta_c, Delta_RF].');
end
end
