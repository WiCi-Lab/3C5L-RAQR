function Omega = calculate_rabi(P, r, mu, Z_0, h_bar)
Intensity = P / (pi * r^2);
Omega = mu * sqrt(2 * Z_0 * Intensity) / h_bar;
end
