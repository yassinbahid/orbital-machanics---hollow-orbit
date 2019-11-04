function V = jacobi_V(x,y,C)
%Earth-Moon: m1=5.98e24;,m2=7.36e22;,mu=m2/(m1 + m2);
global mu
mu = 1/82.45;
V = x.^2+y.^2 +2*(1-mu)./sqrt((x+mu).^2 + y.^2) + 2*mu./sqrt((x+mu-1).^2 + y.^2) - C;