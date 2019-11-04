function J = jacobi_C(X,C)
%Earth-Moon: m1=5.98e24;,m2=7.36e22;,mu=m2/(m1 + m2);
global U
mu = U;
if nargin ==1
    C = 0;
end

x = X(:,1);
y = X(:,2);
vx = X(:,3);
vy = X(:,4);

V = x.^2+y.^2 +2*(1-mu)./sqrt((x+mu).^2 + y.^2) + 2*mu./sqrt((x+mu-1).^2 + y.^2) - C;
J = V - vx.*vx - vy.*vy;