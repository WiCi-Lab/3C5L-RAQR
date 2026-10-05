function [rho_ss_matrix, rho_ss_vec] = solve_steady_state_5L(Op, Od, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec)
L_matrix = zeros(25, 25); eye25 = eye(25);
for k = 1:25
    L_matrix(:, k) = bloch_5L_matrix(0, eye25(:, k), Op, Od, Oc, Olo, Gamma_pop, Gamma_coh, gamma_t, D_vec);
end
trace_eq = zeros(1, 25); trace_eq([1, 7, 13, 19, 25]) = 1;
L_matrix(25, :) = trace_eq;
b = zeros(25, 1); b(25) = 1;
rho_ss_vec = L_matrix \ b;
rho_ss_matrix = reshape(rho_ss_vec, 5, 5);
end
