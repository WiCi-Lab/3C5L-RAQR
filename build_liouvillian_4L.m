function L=build_liouvillian_4L(Op,Oc,Olo,Gamma_pop,Gamma_coh,gamma_t,D_vec)
L=zeros(16,16); E=eye(16);
for k=1:16
    L(:,k)=bloch_4L_matrix(0,E(:,k),Op,Oc,Olo,Gamma_pop,Gamma_coh,gamma_t,D_vec);
end
end
