clear, close all

%% Parameter values
eps = 0.01;             % 1/distance between resonators
N = 4;                  % number of resonators


%%
% initial guess for the N-1 parameters and the single eigenvalue (gamma)
x0 = zeros(1,N);
for n = 1:N
    if n == N
        x0(n) = 1;              % guess 1 for gamma
    elseif mod(n,2) == 0
        x0(n) = 1;              % guess 1 for the real parts
    else
        x0(n) = eps;            % guess epsilon for the imaginary parts
    end
end

f = @(x) excep(N,x(1:end-1),x(end),eps);
options = optimoptions('fsolve','display','iter','MaxFunctionEvaluations',n*1e3);
soln = fsolve(f,x0,options);

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

