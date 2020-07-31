%%% Initializes kappab which to leading order in eps corresponds to an
%%% excpetional point

function kappaEP = makeKappaEP(N)
kappaEP = zeros(1,N);
if N == 3 %% Numbers from asyptotics in Mathematica
    a1s = [0, 0.4832780319391283];
    b1s = [1.5257301701322068, 0];
elseif N == 4
    a1s = [0, 0.653857921541393];
    b1s = [1.8727446815685915, 0.5636519844273622];
    %a1s = [0, -0.862898];
    %b1s = [0.0455948, 1.99533];
elseif N == 5
    a1s = [0, 0.7433686646284927, 0.9108533182462492];
    b1s = [2.1277864103331647, 0.948592103148598, 0];
elseif N == 6
    a1s = [0, 0.7995091404934966, 1.053239878393793];
    b1s = [2.327472034904844, 1.2368578609510557, 0.39457913027490504];
elseif N == 7
    a1s = [0, 0.8385269915559069, 1.1453819914774919, 1.2310017875661599];
    b1s = [2.4904839569124615, 1.465061374173415, 0.6949570191530565, 0];
elseif N == 8
    a1s = [0, 0.8675182070472637, 1.2106740992781007, 1.3486427522363564];
    b1s = [2.6275619089267157, 1.6526158002508826, 0.9354747163943931, 0.3042088404629951];
else
    disp('error')
end
      
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

end