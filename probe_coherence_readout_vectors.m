function [c_re,c_im] = probe_coherence_readout_vectors(n)

N = n^2;
idx21 = 2;
idx12 = 1 + n;

c_re = zeros(1,N);
c_im = zeros(1,N);

c_re(idx21) = 1/2;
c_re(idx12) = 1/2;

c_im(idx21) = 1/(2i);
c_im(idx12) = -1/(2i);
end
