function b = frac_drift(U,dU,x,a)
% Method from the paper <Fractional Langevin Monte Carlo...>
% Stable implementation of drift term to avoid 0/0 error for large x
% INPUT: x = scalar or column vector; a = noise expo; type = well type
% U,dU are function handles

a = a - 2; %Levy exponent converted to that used for Riesz derivative

h = 0.1; %discretization
L = 8; %spatial range
K = floor(L/h+1);
k = -K:K; %row vector
g = (-1).^k.*gamma(a+1)./gamma(a/2 -k+1)./gamma(a/2 +k +1);

l = U(x)-U(x-k*h);
l2 = max(l,[],2);
f = -dU(x-k*h).*exp(l-l2);
%f(f<0)=0;
b = sum(g.*f,2)/(h^a).*exp(l2);

%clipping
b(b>1e3) = 500;%1e3; %limit the height of the well
b(b<-1e3) = -500;%1e3; % this is set such that abs(b*dt) = 1