clear, close all

%% Parameter values
eps = 0.01;             % 1/distance between resonators
N = 5;                  % number of resonators


%%
% initial guess for the N-1 parameters and the single eigenvalue (gamma)
% Note: we actually solve for (real part-1) and (gamma-1) so that the
% variables are all of similar order
x0 = zeros(1,N);
for n = 1:N
    if n == N
        x0(n) = 0;              % guess: gamma = 1
    elseif mod(n,2) == 0
        x0(n) = 0;              % guess: real parts = 1
    else
        x0(n) = 5*(N+1-n)*eps;            % guess: imaginary parts = O(epsilon)
    end
end

re = zeros(1,N-1);
if mod(N,2) == 0
    re(2:2:end-1) = 1;
else
    re(2:2:end) = 1;
end

f = @(x) excep(N,re+x(1:end-1),1+x(end),eps);
options = optimoptions('fsolve','display','iter',...
    'MaxFunEvals',N*1e3,'TolFun',1e-6);
[soln, fval, exitflag] = fsolve(f,x0,options);

soln = soln + [re,1];

disp('-- NUMERICAL VALUES --')
printresults(soln)

%% Asymptotic formulas, for comparison (where possible)

if N == 3
    p = [-8/27,0,-2,1];
    out = roots(p);
    c1 = min(abs(out));
    b1 = sqrt(9/4+c1^2/3);
    asoln = [b1*eps, 1+c1*eps, 1+c1*eps/3];
    disp('-- ASYMPTOTIC VALUES --')
    printresults(asoln)
elseif N == 2
    asoln = [eps/sqrt(1-eps^2),1];
    disp('-- ASYMPTOTIC VALUES --')
    printresults(asoln)
end

