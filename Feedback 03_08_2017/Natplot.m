function U=Natplot(L2,mu),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st step: Finding the linearized equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u=mu;
    FL2 = Jmatrix(L2(1),L2(2),mu);
    [P,D] = eig(FL2);
    a = real(D(1,1)),
    b = imag(D(1,1));
    for i= 1:3
        if a ==0 && b ~= 0
            k = i;
            Vr= real( P(:,i));
            Vi = imag( P(:,i));
            break
        end
        a = real(D(i+1,i+1));
        b = imag(D(i+1,i+1));
    end
    b=abs(b);
    % pick g and d.
%     Vi= -[-0.0000;0.3790; 0.2467;0];
% %     Vr=[ 0.1056;-0.0000;-0.0000; -0.8856];
%     Vi=Vi;
    Vr=-Vr;
    

    g=.001; 
    d=.001;
    S1S =@(t) 2*((g*Vr-d*Vi).*cos(b*t)+(g*Vr+d*Vi).*sin(b*t))+[L2(1),0,0,0]';
    S1 =@(t) 2*((g*Vr-d*Vi).*cos(b*t)+(g*Vr+d*Vi).*sin(b*t));
    S1d =@(t) 2*(-b*(g*Vr-d*Vi).*sin(b*t)+b*(g*Vr+d*Vi).*cos(b*t));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 2ns step: Finding the control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      

    r13=@(x,y) ((x+mu)^2+y^2)^(1/2);
    
    r23 = @(x,y) ((x+mu-1)^2+y^2)^(1/2);    
     Nx =@(x,y,vx,vy) x+2*vy-(1-mu)*(x+mu)/r13(x,y)^3-mu*(x+mu-1)/r23(x,y)^3;
    
    Ny =@(x,y,vx,vy)  y-2*vx-(1-mu)*y/r13(x,y)^3-mu*y/r23(x,y)^3;
    
    N=@(x,y,vx,vy) [0,0,-Nx(x,y,vx,vy),-Ny(x,y,vx,vy)]';
    
    S=@(t,X) [X(3),X(4),0,0]'+ N(X(1),X(2), X(3),X(4))
    [T,S] = ode45(S,[0,100*pi/b], S1S(0));
    
         plot(S(:,1),S(:,2),'b',L2(1),0,'g*');
