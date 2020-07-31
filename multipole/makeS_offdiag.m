function M = makeS_offdiag(ci,cj,Rj,k0,N_multi,Jdata_k0Ri,Jdata_k0Rj)
% make the off-diagonal blocks of S

Sk0=zeros(N_multi+1);
const=-1i*k0*Rj^2;

r_x = norm(cj-ci);

for l=1:N_multi+1
for lp=1:N_multi+1
   Jlpk0= Jdata_k0Rj(lp);
   Jlk0= Jdata_k0Ri(l);   
   Sk0(l,lp)=const*Jlpk0*A_coeff(lp-1,l-1,N_multi,k0*r_x)*Jlk0;   
end
end

M=Sk0;

end