%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B
%
% Computes the resonant frequencies of a pair of bubbles with complex
% material coefficients. 
%
% Details of the method are given in the appendices of Ammari, Davies,
% Hiltunen & Yu (2019) "Topologically protected edge modes in..."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

M = 20;
res_store = zeros(3,M);
tau_vals = linspace(0,0.5,M);
for m = 1:M
tau = tau_vals(m);

N = 3;                      % number of domain components / bubbles

%%% Material parameters
% rho1, kappa1, v1 : bubble 1
% rho2, kappa2, v2 : bubble 2
% rho0, kappa0, v0 : background
high = 5000;
   
rho1 = 1;
rho2 = rho1;
rho3 = rho1;
rho0 = high;

% tau = 0.1;
kappa1 = 1+1i*tau;
kappa2 = 1;
kappa3 = 1-1i*tau;
kappa0 = high;

v1 = sqrt(kappa1/rho1);
v2 = sqrt(kappa2/rho2);
v3 = sqrt(kappa3/rho3);

% High contrast parameters \delta
delta=rho1/rho0;            % delta = 1/high

R = [1, 1, 1];               % vector of radii
d = 10;                      % separation between resonators' centres
cx = [0 R(1)+R(2)+d R(1)+2*R(2)+2*d+R(3)];
cy = zeros(1,N);
cz = zeros(1,N);

% Maximum order for expansion (n = 0, +1, +2, ..., +N_multi)
N_multi = 3;

%% Compute initial guesses for the resonances
% Uses a search algorithm to find initial guesses
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) min(eig((MakeA(R,z,rho0,rho1,rho2,rho3,kappa0,kappa1,kappa2,kappa3,delta,N_multi,cx,cy))));

x = linspace(0.02, 0.03, 500);      
init = [];
for correction = [0 0.001i 0.002i 0.003i]
    y = zeros(1, length(x));
    for i = 1:length(x)
        y(i) = abs(f(x(i) - correction));
    end
    for i = 2:length(x)-1
        if y(i)<y(i-1) & y(i)<y(i+1) & (isempty(init) || min(abs(init-x(i)*ones(1,length(init)))) > 1e-8)
            init = [init, x(i) - correction];
        end
    end
end

if length(init) < length(R)
    disp('WARNING: fewer than 2 initial guesses created')
end

init = sort(init);

%% Use Muller's method to compute the resonances

distTol = 5e-5; fTol = 1e-5; iterMax = 10;
resonances = [];
n = 1;
for initGuess = init
        
    z0 = initGuess;
    z1 = initGuess - 0.00001i;
    z2 = initGuess + 0.00001i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e-6
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi \n'], real(res), imag(res))
       resonances = [resonances res];
       n = n + 1;
    end
end

if length(resonances) == 3
    res_store(:,m) = resonances;
elseif length(resonances) == 2
    res_store(:,m) = [resonances, resonances(end)];
end

end

%%

h = plot(tau_vals/high,real(res_store(1,:)),'b',...
    tau_vals/high,real(res_store(2,:)),'b',...
    tau_vals/high,real(res_store(3,:)),'b',...
    tau_vals/high,imag(res_store(1,:)),'r-.',...
    tau_vals/high,imag(res_store(2,:)),'r-.',...
    tau_vals/high,imag(res_store(3,:)),'r-.');
leg = legend(h([1 4]),'Real part','Imaginary part','interpreter','latex','Location','east');
xlabel('Gain/Loss $b$','interpreter','latex')
ylabel('Frequency $\omega$','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)

set(leg,'Position',get(leg,'Position')+[0,0.1,0,0])

