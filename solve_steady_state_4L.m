function [rho_ss_matrix, rho_ss_vec] = solve_steady_state_4L(Op, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec)
L_matrix = zeros(16, 16); eye16 = eye(16);
for k = 1:16
    L_matrix(:, k) = bloch_4L_matrix(0, eye16(:, k), Op, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec);
end
trace_eq = zeros(1, 16); trace_eq([1, 6, 11, 16]) = 1;
L_matrix(16, :) = trace_eq;
b = zeros(16, 1); b(16) = 1;
rho_ss_vec = L_matrix \ b;
rho_ss_matrix = reshape(rho_ss_vec, 4, 4);
end
