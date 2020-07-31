clear

N = 3;
eps = 0.01;

f = @(x) max(abs(diff(eig(Cdv(N,x,eps)))));
g = @(x) f([1+1i*x(1), 1+x(2)]);
x0 = [0,0];
options = optimoptions('fsolve');
out = fsolve(f,x0,options)
eig(Cdv(N,[1+1i*out(1),1+out(2)],eps))









return
%%

N = 2;
eps = 0.1;
b = eps/sqrt(1-eps^2);

A = Cdv(N,1+1i*b,eps);
[V,D] = eig(A)

return


%%

N = 3;          % number of resonators

eps = 0.01;
b = 3*eps/2;
c = 1;

vdel = [1+1i*b, c];        % material parameters on left side

A = Cdv(N,vdel,eps);

eig(A)

return
%%

eps = 0.1;
b_asym = 3*eps/2;

b_vals = linspace(0.001*b_asym,2*b_asym,50);
y = zeros(3,length(b_vals));
for n = 1:length(b_vals)
    b = b_vals(n);
    vdel = [1+1i*b, 1];
    A = Cdv(3,vdel,eps);
    y(:,n) = -1i*sort(1i*eig(A),'ComparisonMethod','real');
end

for n = 1:3
    subplot(3,1,n)
    plot(b_vals,real(y(n,:)),b_vals,imag(y(n,:)))
end

%% 
clear

N = 4;          % number of resonators

eps = 0.01;
b = 3*eps/2;
c = 1;
d = 2*eps;

vdel = [1+1i*b, c+1i*d];        % material parameters on left side

A = Cdv(N,vdel,eps);

eig(A)
