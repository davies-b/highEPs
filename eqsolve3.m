clear

eps = 0.051;

f = @(b,c) [-1/3*(c+2)^2+1+b^2+2*c-eps^2/4*(1+b^2+8*c);...
    -1/27*(c+2)^3+c*(1+b^2)*(1-9*eps^2/4-eps^3)];

g = @(x) f(x(1),x(2));

x0 = [0,1];
out = fsolve(g,x0)