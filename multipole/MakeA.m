function A = MakeA(R,omega,rho0,rho1,rho2,rho3,kappa0,kappa1,kappa2,kappa3,delta,N_multi,cx,cy)

N = length(cx);

k0=omega*sqrt(rho0/kappa0);
k1=omega*sqrt(rho1/kappa1);
k2=omega*sqrt(rho2/kappa2);
k3=omega*sqrt(rho3/kappa3);

Jdata_k0 = [];
Jdata_k1 = [];
Jdata_k2 = [];
Jdata_k3 = [];
Hdata_k0 = [];
Hdata_k1 = [];
Hdata_k2 = [];
Hdata_k3 = [];
dJdata_k0 = [];
dJdata_k1 = [];
dJdata_k2 = [];
dJdata_k3 = [];
dHdata_k0 = [];

for i = 1:N
Jdata_k0 = [Jdata_k0, makeBesselJdata(N_multi, k0*R(i)  )];
Jdata_k1 = [Jdata_k1, makeBesselJdata(N_multi, k1*R(i)  )];
Jdata_k2 = [Jdata_k2, makeBesselJdata(N_multi, k2*R(i)  )];
Jdata_k3 = [Jdata_k3, makeBesselJdata(N_multi, k3*R(i)  )];
Hdata_k0 = [Hdata_k0, makeHankel1data(N_multi, k0*R(i)  )];
Hdata_k1 = [Hdata_k1, makeHankel1data(N_multi, k1*R(i)  )];
Hdata_k2 = [Hdata_k2, makeHankel1data(N_multi, k2*R(i)  )];
Hdata_k3 = [Hdata_k3, makeHankel1data(N_multi, k3*R(i)  )];
dJdata_k0 = [dJdata_k0, makeDeriBesselJdata(N_multi,k0*R(i),Jdata_k0(:,i))];
dJdata_k1 = [dJdata_k1, makeDeriBesselJdata(N_multi,k1*R(i),Jdata_k1(:,i))];
dJdata_k2 = [dJdata_k2, makeDeriBesselJdata(N_multi,k2*R(i),Jdata_k2(:,i))];
dJdata_k3 = [dJdata_k3, makeDeriBesselJdata(N_multi,k3*R(i),Jdata_k3(:,i))];
dHdata_k0 = [dHdata_k0, makeDeriHankel1data(N_multi,k0*R(i),Hdata_k0(:,i))];
end


%%


SpecialFuncDataM = [];
SpecialFuncDataL = [];
for i = 1:N
SpecialFuncDataM = [SpecialFuncDataM, ...
    Jdata_k0(:,i),Jdata_k1(:,i),Jdata_k2(:,i),Jdata_k3(:,i),...
    Hdata_k0(:,i),Hdata_k1(:,i),Hdata_k2(:,i),Hdata_k3(:,i),...
    dJdata_k1(:,i),dJdata_k2(:,i),dJdata_k3(:,i),...
    dHdata_k0(:,i)];
SpecialFuncDataL =[SpecialFuncDataL, Jdata_k0(:,i),dJdata_k0(:,i)];
end


N_oneblock = (N_multi+1)*2; 

A = zeros(N_oneblock*N);

for i=1:N
    for j=1:N
        if i==j
            A((i-1)*N_oneblock+1:i*N_oneblock,(i-1)*N_oneblock+1:i*N_oneblock)= makeA_diag(R(i),k0,k1,k2,k3,delta,N_multi,SpecialFuncDataM(:,12*(i-1)+1:12*i),i);
        else
            A((i-1)*N_oneblock+1:i*N_oneblock,(j-1)*N_oneblock+1:j*N_oneblock)= makeA_offdiag(i,j,[cx(i),cy(i)],[cx(j),cy(j)],R(j),R(i),k0,delta,N_multi,SpecialFuncDataL(:,2*(i-1)+1:2*i),SpecialFuncDataL(:,2*(j-1)+1:2*j));
        end
    end
end

end