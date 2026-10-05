function P_out = get_P_out_calculated(P_in, chi, lam, d)
alpha = (2 * pi / lam) * max(-imag(chi), 0);
P_out = P_in * exp(-alpha * d);
end
