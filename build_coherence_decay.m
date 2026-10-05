function Gamma_coh = build_coherence_decay(Gamma_pop, Gamma_d)
N = numel(Gamma_pop);
Gamma_coh = zeros(N,N);
for i = 1:N
    for j = 1:N
        Gamma_coh(i,j) = (Gamma_pop(i)+Gamma_pop(j))/2;
        if i ~= j
            Gamma_coh(i,j) = Gamma_coh(i,j)+Gamma_d;
        end
    end
end
end
