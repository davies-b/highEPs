%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B. and Hiltunen, E.O.
%
% Computes the leading order approximation of the resonant frequencies of
% an array of resonators with complex material coefficients. Plots the
% frequencies as the gain/loss increases from 0 across an asymptotic exceptional point
%
% Details of the method are given in Ammari et al. (2020) "High-order
% exceptional points and enhanced sensing in subwavelength resonator arrays"
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

%%% Resonator geometry
N = 2;                       % number of domain components / bubbles
even = mod(N,2) == 0;
R = ones(1,N);               % vector of radii
eps = 0.03;

%%% Table of asymptotic EPs 
if N == 2
    a1s = [0];
    b1s = [1];
elseif N == 3
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
elseif N == 5
    a1s = [0, 0.7433686646284927, 0.9108533182462492];
    b1s = [2.1277864103331647, 0.948592103148598, 0];
elseif N == 6
    a1s = [0, 0.7995091404934966, 1.053239878393793];
    b1s = [2.327472034904844, 1.2368578609510557, 0.39457913027490504];
%     
%     a1s = [0, 1.1261613883078554, 2.2685379181629015];
%     b1s = [1.8426158523313096, -2.0774764216154256, 1.1663480769378682];
elseif N == 7
    a1s = [0, 0.8385269915559069, 1.1453819914774919, 1.2310017875661599];
    b1s = [2.4904839569124615, 1.465061374173415, 0.6949570191530565, 0];
elseif N == 8
    a1s = [0, 0.8675182070472637, 1.2106740992781007, 1.3486427522363564];
    b1s = [2.6275619089267157, 1.6526158002508826, 0.9354747163943931, 0.3042088404629951];

%     a1s = [0, 1.71017, 1.08543, 2.71187];
%     b1s = [2.16684,-1.86197, 1.99683, -1.10729];
%     
%     a1s = [0, -2.337573, -0.524136, -1.220064];
%     b1s = [0.455674, -2.492540, 1.440055,-2.017342];
elseif N == 14
    a1s = [0, 0.956620, 1.396629, 1.654193, 1.813961, 1.908746, 1.953199];
    b1s = [3.163613, 2.354460, 1.797260, 1.337931, 0.928904, 0.547786, 0.181099];
    
%     a1s = [0, -1.004371, -1.299028, -3.281664, -0.020124, -1.739938, -0.484783];
%     b1s = [0.223134, -2.003327, 2.024189, -2.997677, 0.897483, -2.398478, 1.754378];
else
    disp('error: no tabulated EP data')
end

%%% Material parameters
% rhob, kappab : resonators (vector)
% rho0, kappa0 : background
high = 5000;
   
rhob = ones(1,N);
rho0 = high;
delta = rhob(1)/rho0;

kappa0 = high;

% Maximum order for expansion (n = 0, +1, +2, ..., +N_multi)
N_multi = 3;

% Loop parameters
M = 1000; % Assume M is even 
tau_vals = linspace(0,2,M);

%% Compute the resonances
res_store = zeros(N,M);
for m = 1:M
    tau = tau_vals(m);
    vdel_1 = a1s + tau*1i*b1s;
    A = Cd1v(N, vdel_1);
    eval = eig(A);
    resonances = sqrt(3*delta*(1+eps*eval));
    res_store(:,m) = resonances;
end

% Sorting the resonances to make the bands continuous
for m = 2:M/2
    I = 1:N;
    for i = 1:N
        [~,I(i)] = min(res_store(i,m-1)-res_store(:,m));
    end
    res_store(:,m) = res_store(I,m);
end
for m = M-1:-1:M/2
    I = 1:N;
    for i = 1:N
        [~,I(i)] = min(res_store(i,m+1)-res_store(:,m));
    end
    res_store(:,m) = res_store(I,m);
end


%% Plotting
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
print(strcat("N=",num2str(N)),'-depsc')
end
