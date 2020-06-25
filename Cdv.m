function B = Cdv(N, vdel, eps)
% Generates the dilute capacitance matrix $C_d^v$ for N resonators with
% vector of parameters vdel and separation distance 1/eps.

if mod(N,2) == 0
    A = diag([vdel, fliplr(conj(vdel))]);
else
    A = diag([vdel, fliplr(conj(vdel(1:end-1)))]);
end

% A = diag(vdel)

if size(A) ~= [N N]
    disp('[Cdv.m] Warning: incorrect dimensions')
end

B = zeros(N,N);
for n = -N+1:N-1
    if n == 0
        B = B + eye(N);
    else
        B = B + diag(-eps/abs(n)*ones(N-abs(n),1),n);
    end
end

B = A*B;

end