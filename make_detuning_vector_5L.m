function D_vec = make_detuning_vector_5L(v,kp,kd,kc,detuning)
if isscalar(detuning)
    d=[0,0,0,detuning];
elseif numel(detuning)==4
    d=reshape(detuning,1,[]);
else
    error('5L detuning must be scalar or [Delta_p, Delta_d, Delta_c, Delta_RF].');
end
D_vec=[d(1)-kp*v, ...
    d(1)+d(2)-kp*v+kd*v, ...
    d(1)+d(2)+d(3)-kp*v+kd*v-kc*v, ...
    d(1)+d(2)+d(3)+d(4)-kp*v+kd*v-kc*v];
end
