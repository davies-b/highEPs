clear, close all
% Attempt to find the points where the eigenvalues of C^v coalesce numerically

%% Resonator geometry
N = 3;                       % number of domain components / bubbles
R = ones(1,N);               % vector of radii
eps = 0.01;
d = 1/eps;                   % separation between resonators' centres
cx = zeros(1,N);
for i = 2:N
    cx(i) = cx(i-1) + d;
end
cy = zeros(1,N);
cz = zeros(1,N);

% Maximum order for expansion (n = 0, +1, +2, ..., +N_multi)
N_multi = 3;

%%
% initial guess for the N-1 parameters and the single eigenvalue (gamma)
% Note: we actually solve for (real part-1) and (gamma-1) so that the
% variables are all of similar order
if N == 3 %% Numbers from dilute asyptotics in Mathematica
    a1s = [0, 0.4832780319391283];
    b1s = [1.5257301701322068, 0];
elseif N == 4
    a1s = [0, 0.653857921541393];
    b1s = [1.8727446815685915, 0.5636519844273622];
    %a1s = [0, -0.862898];
    %b1s = [0.0455948, 1.99533];
elseif N == 5
    a1s = [0, 0.7433686646284927, 0.9108533182462492];
    b1s = [2.1277864103331647, 0.948592103148598, 0];
elseif N == 6
    a1s = [0, 0.7995091404934966, 1.053239878393793];
    b1s = [2.327472034904844, 1.2368578609510557, 0.39457913027490504];
elseif N == 7
    a1s = [0, 0.8385269915559069, 1.1453819914774919, 1.2310017875661599];
    b1s = [2.4904839569124615, 1.465061374173415, 0.6949570191530565, 0];
elseif N == 8
    a1s = [0, 0.8675182070472637, 1.2106740992781007, 1.3486427522363564];
    b1s = [2.6275619089267157, 1.6526158002508826, 0.9354747163943931, 0.3042088404629951];
else
    disp('error')
end
x0 = zeros(1,N);
for n = 1:N-1
    if mod(n,2) == 0
        x0(n) = a1s(n/2+1)*eps;              % guess: real parts = 1
    else
        x0(n) = b1s((n+1)/2)*eps;            % guess: imaginary parts = O(epsilon)
    end
end
if mod(N,2) == 0
    x0(N) = mean(a1s)*eps;              % guess: gamma = mean of a_i
else
    x0(N) = mean([a1s, a1s(1:end-1)])*eps;
end

re = zeros(1,N-1);
if mod(N,2) == 0
    re(2:2:end-1) = 1;
else
    re(2:2:end) = 1;
end
excep(R,re+x0(1:end-1),1+x0(end),N_multi,cx,cy);

f = @(x) excep(R,re+x(1:end-1),1+x(end),N_multi,cx,cy);
options = optimoptions('fsolve','display','iter',...
    'MaxFunEvals',N*1e3,'TolFun',1e-6);
[soln, fval, exitflag] = fsolve(f,x0,options);
soln_1 = soln/eps;
soln = soln + [re,1];

disp('-- NUMERICAL VALUES --')
printresults(soln)
