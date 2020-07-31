function A = MakeS(R,omega,rho0,kappa0,N_multi,cx,cy)

N = length(cx);

k0=omega*sqrt(rho0/kappa0);

Jdata_k0 = [];
Hdata_k0 = [];

for i = 1:N
    Jdata_k0 = [Jdata_k0, makeBesselJdata(N_multi, k0*R(i)  )];
    Hdata_k0 = [Hdata_k0, makeHankel1data(N_multi, k0*R(i)  )];
end

%%
N_oneblock = N_multi+1; 
A = zeros(N_oneblock*N);

for i=1:N
    for j=1:N
        if i==j
            A((i-1)*N_oneblock+1:i*N_oneblock,(i-1)*N_oneblock+1:i*N_oneblock)= makeS_diag(R(i),k0,N_multi,Jdata_k0(:,i),Hdata_k0(:,i));
        else
            A((i-1)*N_oneblock+1:i*N_oneblock,(j-1)*N_oneblock+1:j*N_oneblock)= makeS_offdiag([cx(i),cy(i)],[cx(j),cy(j)],R(j),k0,N_multi,Jdata_k0(:,i),Jdata_k0(:,j));
        end
    end
end

end