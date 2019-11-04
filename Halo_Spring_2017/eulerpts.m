function L = eulerpts(mu,Jacobi)
% Computes Euler colinear libration points
% by solving the quintic.  
% mu is primaries mass ratio m_2/M = m_2/(m_1+m_2)
% Earth-Moon: m1=5.98e24;,m2=7.36e22;,mu=m2/(m1 + m2);

global MU
MU = mu;

euler_h = @eulerptspoly;
bounds=[-1-MU, -(1+1.e-5)*MU];
L(3) =fzero(euler_h,bounds);


bounds=[0, 1.-1.e-5-MU];
L(1) =fzero(euler_h,bounds);

bounds=[1.0+1.e-5-MU, 2];
L(2) =fzero(euler_h,bounds);

L(4) = .5-MU;
L(5) = -sqrt(3)/2;

if nargin > 1
    wid=1.7;
    low=-2.;
    xj=-wid:.02:wid;
    yj = xj;
    [XJ,YJ]=meshgrid(xj,yj);
    P = -PCR3B_eff_pot(XJ,YJ);
    lows=find(P<low);
    P(lows)=low;
    figure
    C=contour(XJ, YJ, -P,30);
    clabel(C,'manual')
    hold
    ptlab = 'bp';
    plot(L(1),0,ptlab)
    ptlab = 'mp';
    plot(L(2),0,ptlab)
    ptlab = 'cp';
    plot(L(3),0,ptlab)
    ptlab = 'r*';
    plot(L(4),L(5),ptlab)
    ptlab = 'r*';
    plot(L(4),-L(5),ptlab)
    plot(-MU,0,'bo')
    plot(1-MU,0,'yo')
    xlabel('x')
    ylabel('y')
    title('Equipotentials of PCR3B:   \Omega = J - K <= J')
    
    figure
    surf(XJ,YJ,P)
    title('PCR3B Effective Potential Surface z=-\Omega')
    xlabel('x')
    ylabel('y')
    zlabel('-\Omega (Effective Potential)')
    hold on
    top1=-1.5;
    plot3([L(1),L(1)],[0,0],[low,top1],'k-',L(1),0,top1,'bp')
    top2=-1.5;
    plot3([L(2),L(2)],[0,0],[low,top2],'k-',L(2),0,top2,'mp')
    top3=-1.45;
    plot3([L(3),L(3)],[0,0],[low,top3],'k-',L(3),0,top3,'cp')
    top4=-1.45;
    plot3([L(4),L(4)],[-L(5),-L(5)],[low,top4],'k-',L(4),-L(5),top4,'r*')
    plot3([L(4),L(4)],[ L(5), L(5)],[low,top4],'k-',L(4), L(5),top4,'r*')
end

end

function P = PCR3B_eff_pot(x,y)
% Computes effective potential P 
% at location(s) (x,y) 
% for Planar Circular Restricted 3-Body Problem
% with primaries mass ratio MU a global variable

global MU


r1 = sqrt((x+MU).^2 + y.^2);
r2 = sqrt((x -1+MU).^2 + y.^2);
P = .5*(x.*x + y.*y) +(1-MU)./r1 + MU./r2;
end

function y = eulerptspoly(x)
% Evaluates polynomial for Euler's
% Colinear Libration points
% (for use with root finder)

global MU


s1=x+MU;
ss1 = s1.*s1;
sg1 = sign(s1);
sg1(find(sg1==0)) = 1;

s2=x-1+MU;
ss2 = s2.*s2;
sg2 = sign(s2);
sg2(find(sg2==0)) = 1;



y = x.*ss1.*ss2 - (1-MU).*sg1.*ss2 - MU*sg2.*ss1;

end