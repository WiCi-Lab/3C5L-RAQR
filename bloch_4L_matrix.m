function drho_vec = bloch_4L_matrix(~, rho_vec, Op, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec)
rho = reshape(rho_vec, 4, 4); H = zeros(4,4);
H(1,2) = Op/2;  H(2,1) = conj(Op)/2;
H(2,3) = Oc/2;  H(3,2) = conj(Oc)/2;
H(3,4) = Olo/2; H(4,3) = conj(Olo)/2;
H(2,2) = -D_vec(1); H(3,3) = -D_vec(2); H(4,4) = -D_vec(3);

comm = -1i * (H * rho - rho * H);
L_diss = zeros(4,4);
L_diss(4,4) = -(Gamma_pop(4) + gamma_t) * rho(4,4);
L_diss(3,3) =  Gamma_pop(4) * rho(4,4) - (Gamma_pop(3) + gamma_t) * rho(3,3);
L_diss(2,2) =  Gamma_pop(3) * rho(3,3) - (Gamma_pop(2) + gamma_t) * rho(2,2);
L_diss(1,1) =  Gamma_pop(2) * rho(2,2) + gamma_t * (rho(2,2) + rho(3,3) + rho(4,4));
for i = 1:4
    for j = 1:4
        if i ~= j, L_diss(i,j) = -(Gamma_coh(i,j) + gamma_t) * rho(i,j); end
    end
end
drho_vec = reshape(comm + L_diss, [], 1);
end
