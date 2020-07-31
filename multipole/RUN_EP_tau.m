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

clear all, %close all

%%% Resonator geometry
N = 4;                       % number of domain components / bubbles
even = mod(N,2) == 0;
R = ones(1,N);               % vector of radii
eps = 0.1;
%eps = 1/12;
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
if N == 3 %% Numbers from asyptotics in Mathematica
    a1s = [0, 0.4832780319391283];
    b1s = [1.5257301701322068, 0];
elseif N == 4
    a1s = [0, 0.653857921541393];
    b1s = [1.8727446815685915, 0.5636519844273622];
end
    %a1s = [0, 0.467507115978696]; % From optimization of C^v
    %b1s = [1.519904079902665, 0];
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
tau = tau_vals(1);

kappab = real(kappaEP) + 1i*tau*imag(kappaEP);

%% Compute initial guesses for the resonances
% Uses a search algorithm to find initial guesses
% Define function f : f gives minimum of eigenvalues of the operator A
% MakeA : gives a matrix approximation for the operator A

f= @(z) eigs(MakeA(R,z,rho0,rhob,kappa0,kappab,N_multi,cx,cy),1,'smallestabs');
x = linspace(0.022, 0.027, 1000);      
init = [];
for correction = [0 0.00001i 0.00005i 0.0001i 0.0002i 0.00051i]
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

%% Use capacitance to compute initial guesses
% %Cv = 4*pi*MakeCv(R,kappab(1:ceil(N/2)),N_multi,cx,cy);
% Cv = 4*pi*Cdv(N, kappab(1:ceil(N/2)), eps);
% vol = 4*pi*R(1)^3/3;
% init = sort(sqrt(1/high*eig(Cv)/vol)); 
% %init = init.';
% %init = [init, init + 1/high*(1-1i), init.*linspace(1.01,0.99,N)];
% %init = linspace(init(1),init(end),100);

%% Use Muller's method to compute the resonances
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
%leg = legend('Real part','Imaginary part','interpreter','latex','Location','east');
ylabel('Real part','interpreter','latex')
ax = gca;
ax.YAxis.Exponent = -2;
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)
subplot(2,1,2)
%leg = legend('Real part','Imaginary part','interpreter','latex','Location','east');
xlabel('Gain/Loss $\tau$','interpreter','latex')
ylabel('Imaginary part','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)


prt = 1;
if prt 
print(strcat("N=",num2str(N),"full"),'-depsc')
end

return

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