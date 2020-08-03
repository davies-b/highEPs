function out = excep(N, pars, gamma)
%%% Input
% N: number of resonators
% pars: the material parameters (denoted b,c,d,...) in the paper
% gamma: the single eigenvalue
%%% Output
% The difference between then polynomial coefficients of the characteristic
% polynomial of $C_{d,1}^v$ and its desired form: (x-gamma_1)^N

if length(pars)~=N-1
    disp('[excep2.m] Warning: incorrect number of parameters given')
end

vdel_1 = [1i*pars(1)];
n = 2;
while n <= length(pars)
    if n<length(pars)
        vdel_1 = [vdel_1, pars(n)+1i*pars(n+1)];
        n = n + 2;
    elseif n == length(pars)
        vdel_1 = [vdel_1, pars(n)];
        n = n + 1;
    end
end

M = Cd1v(N,vdel_1);

chpoly = charpoly(M);
chpoly = chpoly(2:end);            % highest order coefficient is always 1

goal = zeros(1,N);
for n = 2:N+1
    goal(n-1) = nchoosek(N,n-1)*(-gamma)^(n-1);
end

out = goal-chpoly;