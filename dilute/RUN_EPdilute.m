% Attempt to find the points where the eigenvalues of C_{d,1}^v coalesce numerically
clear, close all

%% Parameter values
N = 4;                  % number of resonators

%%
% initial guess for the N-1 parameters and the single eigenvalue (gamma1)
x0 = zeros(1,N);
Nb = ceil(N/2);
altsign = (-1).^(0:Nb-1);

b1s = altsign; % For example, look for alternating solutions
a1s = zeros(1,Nb);
for n = 1:N-1
    if mod(n,2) == 0
        x0(n) = a1s(n/2+1);              % Inital guesses of a
    else
        x0(n) = b1s((n+1)/2);            % Initial guesses of b
    end
end
if mod(N,2) == 0
    x0(N) = mean(a1s);              % guess: gamma = mean of a_i
else
    x0(N) = mean([a1s, a1s(1:end-1)]);
end

f = @(x) excep(N,x(1:end-1),x(end));
options = optimoptions('fsolve','display','iter',...
    'MaxFunEvals',N*1e3,'TolFun',1e-6);
[soln, fval, exitflag] = fsolve(f,x0,options);

disp('-- NUMERICAL VALUES --')
printresults(soln)

%% Asymptotic formulas, for comparison (where possible)

if N == 3
    p = [-8/27,0,-2,1];
    out = roots(p);
    c1 = min(abs(out));
    b1 = sqrt(9/4+c1^2/3);
    asoln = [b1, c1, c1/3];
    disp('-- ASYMPTOTIC VALUES --')
    printresults(asoln)
elseif N == 2
    asoln = [1,1];
    disp('-- ASYMPTOTIC VALUES --')
    printresults(asoln)
end