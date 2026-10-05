function Hf = evaluate_baseband_model(model,f_Hz)

f_shape = size(f_Hz);
f = f_Hz(:).';
N = model.levels^2;
I = eye(N);
H = zeros(size(f));

for iv = 1:model.Nv
    L0 = model.L0{iv};
    rhs = model.rhs{iv};
    w = model.weights(iv);
    for ik = 1:numel(f)
        A = 1i*2*pi*f(ik)*I - L0;
        A(N,:) = model.trace_eq;
        drho = A\rhs;
        H(ik) = H(ik) + model.K*(model.c_det*drho)*w;
    end
end

Hf = reshape(H,f_shape);
end
