function out = excep(R,pars,gamma,N_multi,cx,cy)
%%% Input
% N: number of resonators
% pars: the material parameters (denoted b,c,d,...) in the paper
% gamma: the single eigenvalue
% eps: 1/distance between the resonaots
%%% Output
% The difference between then polynomial coefficients of the characteristic
% polynomial of $C_d^v$ and its desired form: (x-gamma)^N

N = length(cx);

if length(pars)~=N-1
    disp('[excep.m] Warning: incorrect number of parameters given')
end

vdel = [1+1i*pars(1)];
n = 2;
while n <= length(pars)
    if n<length(pars)
        vdel = [vdel, pars(n)+1i*pars(n+1)];
        n = n + 2;
    elseif n == length(pars)
        vdel = [vdel, pars(n)];
        n = n + 1;
    end
end
A = MakeCv(R,vdel,N_multi,cx,cy);

chpoly = charpoly(A);
chpoly = chpoly(2:end);            % highest order coefficient is always 1

goal = zeros(1,N);
for n = 2:N+1
    goal(n-1) = nchoosek(N,n-1)*(-gamma)^(n-1);
end

out = goal-chpoly;