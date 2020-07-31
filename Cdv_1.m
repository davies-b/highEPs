function M = Cdv_1(N, vdel_1)
% Generates the first order correction to dilute capacitance matrix, in other words $1/\epsilon(C_d^v - \gamma I)$
% for N resonators with vector of parameters vdel and separation distance 1/eps.

if mod(N,2) == 0
    d = [vdel_1, fliplr(conj(vdel_1))];
else
    d = [vdel_1, fliplr(conj(vdel_1(1:end-1)))];
end
M = diag(d);

if size(M) ~= [N N]
    disp('[Cdveps.m] Warning: incorrect dimensions')
end

for i = 1:N
    for j = 1:N
        if i ~= j
            M(i,j) = -1/abs(i-j);
        end
    end
end

end