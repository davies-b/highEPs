function Cv = MakeCv(R,vdel,N_multi,cx,cy)
N = length(cx);
if mod(N,2) == 0
    A = diag([vdel, fliplr(conj(vdel))]);
else
    A = diag([vdel, fliplr(conj(vdel(1:end-1)))]);
end

omega = 0.0001;
M = N_multi+1;
S = MakeS(R,omega,1,1,N_multi,cx,cy);
C = zeros(N,N);
for i = 1:N
    phii = zeros(N*M,1); phii(M*(i-1) + 1) = 1;
    psii = S\phii;
    for j = 1:N
        phij = zeros(N*M,1); phij(M*(j-1) + 1) = 1;
        C(j,i) = -4*pi*R(i)^2*(phij'*psii);
    end
end
C = C/(4*pi); % rescale to make gamma = 1 + O(eps) (leading order 1)
Cv = A*C;
end
