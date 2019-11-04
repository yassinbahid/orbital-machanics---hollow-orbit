
function L2=eulervect(mu,n)
% Program by Bruce Lundberg, PHD, and Yassin Bahid
%computing the second libration vector:
%n is the vectro you want: L1, L2, or L3
% mu is primaries mass ratio m_2/M = m_2/(m_1+m_2)
% Earth-Moon: m1=5.98e24;,m2=7.36e22;,mu=m2/(m1 + m2);
    bounds1 = [1.0+1.e-5-mu, 2];
    bounds2=[0, 1.-1.e-5-mu];
    bounds3=[-1-mu, -(1+1.e-5)*mu];
    Equil = @(x) x-(1-mu)*(x+mu)/(sqrt((x+mu)^2))^3-mu*(x+mu-1)/(sqrt((x+mu-1)^2))^3;
    P1 = fzero(Equil,bounds1);
    P2 = fzero(Equil,bounds2);
    P3 = fzero(Equil,bounds3);
    L1=[P1,0,0,0];
    L2=[P2,0,0,0];
    L3=[P3,0,0,0];
    L=[L1,L2,L3];
    
end
