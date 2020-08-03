%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B., Hiltunen, E.O.
%
% Computes the resonant frequencies of an array of resonators with complex
% material coefficients using the multipole method. Plots the frequencies
% as the gain/loss increases from 0 across an asymptotic exceptional point
%
% Details of the method are given in the appendices of Ammari, Davies,
% Hiltunen & Yu (2019) "Topologically protected edge modes in..."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

%%% Resonator geometry
N = 3;                       % number of domain components / bubbles
even = mod(N,2) == 0;
R = ones(1,N);               % vector of radii
eps = 0.1;
d = 1/eps;                   % separation between resonators' centres
cx = zeros(1,N);
for i = 2:N
    cx(i) = cx(i-1) + d;
end
cy = zeros(1,N);
cz = zeros(1,N);

%%% Material parameters
% rhob, kappab : resonators (vector)
% rho0, kappa0 : background
high = 5000;
   
rhob = ones(1,N);
rho0 = high;

kappa0 = high;
kappaEP = zeros(1,N);
if N == 3 
    a1s = [0, 0.4832780319391283];
    b1s = [1.5257301701322068, 0];
elseif N == 4
    a1s = [0, 0.653857921541393];
    b1s = [1.8727446815685915, 0.5636519844273622];
   
%     a1s = [0, -0.8628981032057736];
%     b1s = [0.045594766758086464, 1.9953267393119942];
%     
%     a1s = [0, 1.0660466602568839];
%     b1s = [1.7017496802269811, -1.132866663769782];
%     
%     a1s = [0, -1.1541795039569698];
%     b1s = [0.7341129919407204, -1.9334565911476076];
else 
    disp('error: no tabulated EP data')
end

for i = 1:N/2
    kappaEP(i) = 1 + a1s(i)*eps + 1i*b1s(i)*eps;
end
if even
    kappaEP(N/2+1:end) = fliplr(conj(kappaEP(1:N/2)));
else
    i = (N+1)/2;
    kappaEP(i) = 1 + a1s(i)*eps;
    kappaEP(i+1:end) = fliplr(conj(kappaEP(1:i-1)));
end

% Maximum order for expansion (n = 0, +1, +2, ..., +N_multi)
N_multi = 3;

% Loop parameters
M = 100;
tau_vals = linspace(0,2,M);

%% Compute initial guesses for the resonances
% Uses a search algorithm to find initial guesses
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

tau = tau_vals(1);
kappab = real(kappaEP) + 1i*tau*imag(kappaEP);

f= @(z) eigs(MakeA(R,z,rho0,rhob,kappa0,kappab,N_multi,cx,cy),1,'smallestabs');
x = linspace(0.022, 0.027, 100);      
init = [];
for correction = [0 0.00001i 0.0001i 0.00051i]
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
    disp('ERROR: fewer than N initial guesses created')
    return
end

init = sort(init);
init0 = init;

%% Use Muller's method to compute the resonances
% Uses the previous values for kappab as initial guesses for the next
% values
res_store = zeros(N,M);
distTol = 5e-5; fTol = 1e-5; iterMax = 20;
for m = 1:M
tau = tau_vals(m);
kappab = real(kappaEP) + 1i*tau*imag(kappaEP);
f= @(z) eigs(MakeA(R,z,rho0,rhob,kappa0,kappab,N_multi,cx,cy),1,'smallestabs');
resonances = [];
n = 1;
for initGuess = init
        
    z2 = initGuess;
    z1 = initGuess - 0.00001i;
    z0 = initGuess + 0.00001i;
    
    res = MullersMethod(f, z0, z1, z2, iterMax, distTol, fTol);
    if isempty(resonances) || min(abs(resonances-res*ones(1,length(resonances)))) > 1e-6
       fprintf(['Resonant frequency #', num2str(n), ' :   %.8f %.8fi \n'], real(res), imag(res))
       resonances = [resonances res];
       n = n + 1;
    end
end

res_store(:,m) = resonances;
init = resonances;
end

%%
figure
for i = 1:N
    subplot(2,1,1)
    hold on
    plot(tau_vals,real(res_store(i,:)),'b');
    subplot(2,1,2)
    hold on
    plot(tau_vals,imag(res_store(i,:)),'r');
end
subplot(2,1,1)
ylabel('Real part','interpreter','latex')
ax = gca;
ax.YAxis.Exponent = -2;
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)
subplot(2,1,2)
xlabel('Gain/Loss $\tau$','interpreter','latex')
ylabel('Imaginary part','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)

% create smaller axes in top right, and plot on it
if even
    b1se = [b1s, -fliplr(b1s)];
else
    b1se = [b1s, -fliplr(b1s(1:end-1))];
end
axes('Position',[.74 .84 .18 .14])
box on
clr = [0.3 0.3 0.3];
bar(b1se,'EdgeColor', clr,'FaceColor', clr)
set(gca,'xtick',[])

prt = 0;
if prt 
    print(strcat("N=",num2str(N),"full"),'-depsc')
end