function printresults_1(soln)

n = 1;
m = 1;
while n < length(soln)
    if n == 1
        fprintf(['Resonantor #', num2str(m), ' :   %.6f + %.6fi \n'], 0, soln(n))
        n = n + 1;
    elseif n < length(soln)-1
        fprintf(['Resonantor #', num2str(m), ' :   %.6f + %.6fi \n'], soln(n), soln(n+1))
        n = n + 2;
    elseif n == length(soln)-1
        fprintf(['Resonantor #', num2str(m), ' :   %.6f \n'], soln(n))
        n = n + 1;
    end
    m = m + 1;
end
fprintf(['Eigenvalue    :   %.6f \n'], soln(end))