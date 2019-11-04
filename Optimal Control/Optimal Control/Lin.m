function [b,S10] = Lin(L2,mu)
% Program by Bruce Lundberg, PHD, and Yassin Bahid
% Finding the linearized equation

    u=mu;
    FL2 = Jmatrix(L2(1),L2(2),mu),
    [P,D] = eig(FL2),
    a = real(D(1,1));
    b = imag(D(1,1));
    for i= 1:3
        if abs(a) < 10*eps && b ~= 0
            k = i;
            Vr= real( P(:,i)),
            Vi = imag( P(:,i));
            break
        end
        a = real(D(i+1,i+1));
        b = imag(D(i+1,i+1));
    end
    b=abs(b);
    % pick g and d.
    g=.001; 
    d=.001;
    S1S =@(t) 2*((g*Vr-d*Vi).*cos(b*t)+(g*Vr+d*Vi).*sin(b*t))+[L2(1),0,0,0]';
    S1 =@(t) 2*((g*Vr-d*Vi).*cos(b*t)+(g*Vr+d*Vi).*sin(b*t));
    S10 = S1S(0);