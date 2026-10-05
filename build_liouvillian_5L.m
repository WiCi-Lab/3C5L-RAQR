function L=build_liouvillian_5L(Op,Od,Oc,Olo,Gamma_pop,Gamma_coh,gamma_t,D_vec)
L=zeros(25,25); E=eye(25);
for k=1:25
    L(:,k)=bloch_5L_matrix(0,E(:,k),Op,Od,Oc,Olo,Gamma_pop,Gamma_coh,gamma_t,D_vec);
end
end
