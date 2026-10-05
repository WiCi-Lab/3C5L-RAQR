function dL = rf_drive_liouvillian_derivative_4L()

n = 4;
dH = zeros(n,n);
dH(3,4) = 1/2;
dH(4,3) = 1/2;
dL = -1i * (kron(eye(n),dH) - kron(dH.',eye(n)));
end
