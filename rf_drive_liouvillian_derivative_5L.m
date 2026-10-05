function dL = rf_drive_liouvillian_derivative_5L()

n = 5;
dH = zeros(n,n);
dH(4,5) = 1/2;
dH(5,4) = 1/2;
dL = -1i * (kron(eye(n),dH) - kron(dH.',eye(n)));
end
