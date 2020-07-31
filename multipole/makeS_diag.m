function M = makeS_diag(R,k0,N_multi,Jdata_k0R,Hdata_k0R)
% make the diagonal blocks of A

const = -1i*R^2;

Sk0=zeros(N_multi+1);

for n=1:N_multi+1
   
   Jk0 = Jdata_k0R(n);
   Hk0 = Hdata_k0R(n);
    
   Sk0(n,n) = const*k0*Jk0*Hk0;
   
end

M=Sk0;
