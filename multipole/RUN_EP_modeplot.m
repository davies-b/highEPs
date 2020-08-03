%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Davies, B. and Hiltunen, E.O.
%
% Plots the eigenmode (squared) for the exceptional points in arrays of 
% three or four resonators, based on the asymptotic values.
%
% Details of the method are given in the appendices of Ammari, Davies,
% Hiltunen & Yu (2020) "Topologically protected edge modes in..."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all

%%% Resonator geometry
N = 4;                       % number of domain components / bubbles
option = 2;                  % picks which EP in the case of N=4
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
kappab = zeros(1,N);
%%% Precomputed EP data
if N == 3
    a1s = [0, 0.4832780319391283];
    b1s = [1.5257301701322068, 0];
    omega = sqrt((3+a1s(2)*eps)/high);
elseif N == 4
    if option == 1
        a1s = [0, 0.653857921541393];
        b1s = [1.8727446815685915, 0.5636519844273622];
    elseif option == 2
        a1s = [0, -0.8628981032057736];
        b1s = [0.045594766758086464, 1.9953267393119942];
    elseif option == 3
        a1s = [0, 1.0660466602568839];
        b1s = [1.7017496802269811, -1.132866663769782];
    elseif option == 4
        a1s = [0, -1.1541795039569698];
        b1s = [0.7341129919407204, -1.9334565911476076];
    end
    omega = sqrt(3/2*(a1s(2)*eps+2)/high);
end

for i = 1:N/2
    kappab(i) = 1 + a1s(i)*eps + 1i*b1s(i)*eps;
end
if even
    kappab(N/2+1:end) = fliplr(conj(kappab(1:N/2)));
else
    i = (N+1)/2;
    kappab(i) = 1 + a1s(i)*eps;
    kappab(i+1:end) = fliplr(conj(kappab(1:i-1)));
end

%% Compute the eigenmode in a plane

N_multi = 0;
A = MakeA(R,omega,rho0,rhob,kappa0,kappab,N_multi,cx,cy);
[V, D] = eig(A);
[~,permutation] = sort(diag(D));
D = D(permutation,permutation); 
V = V(:,permutation);

phi = [];
psi = [];
N_terms = 1;

for i = 1:N
    phi = [phi, V((i-1)*2*N_terms+1:i*2*N_terms-N_terms,1)];
    psi = [psi, V((i-1)*2*N_terms+N_terms+1:i*2*N_terms,1)];
end

% Calculate field
green = @(k,x,y,z) -exp(1i*k*sqrt(x.^2+y.^2+z.^2))/4/pi./sqrt(x.^2+y.^2+z.^2);

% Grid for field
gridN = N*100+1;  
gridMinX1 = cx(1)-10-R(1);
gridMaxX1 = cx(end)+R(end)+10;
gridMinX2 = -10;
gridMaxX2 = -gridMinX2;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = linspace(gridMinX2, gridMaxX2, gridN);
[ g1, g2 ] = meshgrid(g1, g2);
gridPoints = [g1(:) g2(:)]';
gridPointsN = length(gridPoints);
gridPoints = [gridPoints; zeros(1,gridPointsN)];

u = zeros(gridPointsN, 1);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);
    
    % Determine whether we are inside or outside the domains
    I = (gridPoint(1)*ones(1,length(cx))-cx).^2 + (gridPoint(2)*ones(1,length(cy))).^2 + (gridPoint(3)*ones(1,length(cz))).^2  <= R.^2 ;
    if sum( I ) > 0
        S = 0;
        I = find(I);
        kb = omega*sqrt(rhob(I)/kappab(I));
        fun1 = @(s,t) green(kb, gridPoint(1)-cx(I)-R(I).*sin(s).*cos(t), gridPoint(2)-cy(I)-R(I).*sin(s).*sin(t), gridPoint(3)-cz(I)-R(I).*cos(s)).*phi(I).*harmonicY(0,0,s,t);
        S = integral2(fun1, 0, pi, -pi, pi);
        u(j) = S;
    else
        S = 0;
        k = omega*sqrt(rho0/kappa0);
        for i = 1:N
            fun2 = @(s,t) green(k, gridPoint(1)-cx(i)-R(i).*sin(s).*cos(t), gridPoint(2)-cy(i)-R(i).*sin(s).*sin(t), gridPoint(3)-cz(i)-R(i).*cos(s)).*psi(i).*harmonicY(0,0,s,t);
            S = S + integral2(fun2, 0, pi, -pi, pi);
        end
        u(j) = S;
    end
end
%% Plotting
fsz = 16;

uTotal = reshape(u, [gridN gridN]);
hFig = figure;
if N == 3
    set(hFig, 'Position', [100 100 800 400]);
elseif N == 4
    set(hFig, 'Position', [100 100 800 400]);
end
surf(g1, g2, abs(uTotal.^2), 'edgecolor', 'none'); 
xlabel('$x_1$','interpreter','latex','FontSize', fsz);
ylabel('$x_2$','interpreter','latex','FontSize', fsz,'Rotation',0);
axis([gridMinX1-0.1, gridMaxX1+0.1, gridMinX2-0.1, gridMaxX2+0.1]);
view(0,90);
colormap(flipud(bone))
set(gca,'TickLabelInterpreter','latex','FontSize', fsz)
box on

% create smaller axes in top right, and plot on it
if N == 4
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
end
