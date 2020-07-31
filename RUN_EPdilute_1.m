clear, close all
% Attempt to find the points where the eigenvalues of C_{d,1}^v coalesce numerically

%% Parameter values
N = 14;                  % number of resonators

%%
% initial guess for the N-1 parameters and the single eigenvalue (gamma)
% Note: we actually solve for (real part-1) and (gamma-1) so that the
% variables are all of similar order
x0 = zeros(1,N);
Nb = ceil(N/2);
altsign = (-1).^(0:Nb-1);

a1s = [0, 0.8675182070472637, 1.2106740992781007, 1.3486427522363564]; a1s = [a1s, 1.2*a1s(end)];
b1s = [2.6275619089267157, 1.6526158002508826, 0.9354747163943931, 0.3042088404629951, 0];
b1s = [linspace(3.2,0.2,Nb-1), 0];
a1s = [0, linspace(0.9,2,Nb-1)];
b1s = altsign;
a1s = [0,zeros(1,Nb-1)]; %linspace(2,1,Nb)];
% a1s = [0, 2, 1.71017, 1.08543, 2.71187];
% b1s = [-2, 2.16684,-1.86197, 1.99683, -1.10729];
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

f = @(x) excep_1(N,x(1:end-1),x(end));
options = optimoptions('fsolve','display','iter',...
    'MaxFunEvals',N*1e3,'TolFun',1e-6);
[soln, fval, exitflag] = fsolve(f,x0,options);

disp('-- NUMERICAL VALUES --')
printresults_1(soln)

%% Asymptotic formulas, for comparison (where possible)

if N == 3
    p = [-8/27,0,-2,1];
    out = roots(p);
    c1 = min(abs(out));
    b1 = sqrt(9/4+c1^2/3);
    asoln = [b1, c1, c1/3];
    disp('-- ASYMPTOTIC VALUES --')
    printresults_1(asoln)
elseif N == 2
    asoln = [1,1];
    disp('-- ASYMPTOTIC VALUES --')
    printresults_1(asoln)
end