clear, close all

% Plots how the eigenvalues of the matrix Cdv vary as the gain/loss is
% varied in an array of 3 resonators

eps = 0.01;
N = 3;
c = 1+0.483*eps;

b_vals = eps*linspace(0.01,3,50);
eigs_re = zeros(N,length(b_vals));
eigs_im = zeros(N,length(b_vals));

for n = 1:length(b_vals)
    b = b_vals(n);
    C = Cdv(N,[1+1i*b,c],eps);
    
%     C = [1+1i*b, -eps*(1+1i*b), -eps*(1+1i*b)/2;...
%         -c*eps, c, -c*eps;...
%         -eps*(1-1i*b)/2, -eps*(1-1i*b), 1-1i*b];
    
%     C = [1+1i*b, -eps*(1+1i*b); -eps*(1-1i*b), 1-1i*b];
    
    eigs = eig(C);
    eigs_re(:,n) = sort(real(eigs));
    eigs_im(:,n) = sort(imag(eigs));
    
end

%%
subplot(2,1,1)
plot(b_vals,eigs_re,'b')

ylabel('Real part','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)

subplot(2,1,2)
plot(b_vals,eigs_im,'r')

xlabel('Gain/Loss $b$','interpreter','latex')
ylabel('Imaginary part','interpreter','latex')
set(gca,'ticklabelinterpreter','latex')
set(gca, 'FontSize',16)
