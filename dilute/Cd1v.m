function M = Cd1v(N, vdel_1)
% Generates the first order correction to dilute capacitance matrix, in other words $C_{d,1}^v$
% for N resonators with vector of parameters vdel and separation distance 1/eps.

if mod(N,2) == 0
    d = [vdel_1, fliplr(conj(vdel_1))];
else
    d = [vdel_1, fliplr(conj(vdel_1(1:end-1)))];
end
M = diag(d);

if size(M) ~= [N N]
    disp('[Cdv_1.m] Warning: incorrect dimensions')
end

for i = 1:N
    for j = 1:N
        if i ~= j
            M(i,j) = -1/abs(i-j);
        end
    end
end

end