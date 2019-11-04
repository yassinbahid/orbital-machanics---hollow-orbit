function FL2 = FeedbackL2(L2,mu)
% Earth-Moon: m1=5.98e24;,m2=7.36e22;,mu=m2/(m1 + m2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st step: Finding the linearized equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u=mu;
    FL2 = Jmatrix(L2(1),L2(2),mu);
    [P,D] = eig(FL2);
    a = real(D(1,1)),
    b = imag(D(1,1));
    for i= 1:3
        if a==0 && b ~= 0
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
%     x=L2(1);
%     y=L2(2);
%     mus=(1-u);
%     r1s=((x+u).^2+y.^2);
%     r2s=((x-1+u).^2+y.^2);
%     F_3Lx = 1 - mus*(y.^2 -2*(x + u).^2)./(r1s.^2.5)  -  u*(y.^2 -2*(x + u-1).^2)./(r2s.^2.5);
% 
%     F_3Ly =    3*mus*((x + u).*y)./(r1s.^2.5)  +  3*u*((x + u-1).*y)./(r2s.^2.5);
% 
%     F_4Lx =    3*mus*((x + u).*y)./(r1s.^2.5)  +  3*u*((x + u-1).*y)./(r2s.^2.5);
% 
%     F_4Ly = 1 - mus*((x + u).^2 - 2*y.^2 )./(r1s.^2.5)  -  u*((x + u-1).^2-2*y.^2)./(r2s.^2.5);
%     
    C = 2*(g*Vr-d*Vi); 
    c1 = C(1,1); c2 = C(2,1); c3 = C(3,1); c4 = C(4,1);
    
    S = 2*(g*Vr+d*Vi); s1 = S(1,1); s2 = S(2,1); s3 = S(3,1); s4 = S(4,1);
%     
%     Linx =@(t) -(c1* F_3Lx+c2* F_3Ly+2*c4).*cos(b*t)+ (s1* F_3Lx+s2* F_3Ly+2*s4).*sin(b*t);
%     
%     Liny =@(t) -(c1* F_4Lx+c2* F_4Ly-2*c3).*cos(b*t)+ (s1* F_4Lx+s2* F_4Ly-2*s3).*sin(b*t);
    
    Lin2 =@(t) FL2*S1(t);
    Lin2(.5*pi/b);
    S1d(.5*pi/b);
    r13=@(x,y) ((x+mu)^2+y^2)^(1/2);
    
    r23 = @(x,y) ((x+mu-1)^2+y^2)^(1/2);
    
    Nx =@(x,y,vx,vy) x+2*vy-(1-mu)*(x+mu)/r13(x,y)^3-mu*(x+mu-1)/r23(x,y)^3;
    
    Ny =@(x,y,vx,vy)  y-2*vx-(1-mu)*y/r13(x,y)^3-mu*y/r23(x,y)^3;
    N=@(x,y,vx,vy) [0,0,-Nx(x,y,vx,vy),-Ny(x,y,vx,vy)]';
%   
%     Ux =@(t,x,y,vx,vy) Linx(t)-Nx(x,y,vx,vy);
%     Uy =@(t,x,y,vx,vy) Liny(t)-Ny(x,y,vx,vy);

    D=diag([0,0,1,1]);
    S1dxy=@(t)(D*S1d(t));
    
    Uc=@(t,x,y,vx,vy) S1dxy(t)-N(x,y,vx,vy);
    XP=@(t,X) [X(3),X(4),0,0]'+ N(X(1),X(2), X(3),X(4))+Uc(t,X(1),X(2),X(3),X(4));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%natural form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     XP=@(t,X) [X(3),X(4),0,0]'+ N(X(1),X(2), X(3),X(4));
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %trial 3
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     XP =@(t,X)[X(3),...
%          X(4),...
%          Nx(X(1),X(2),X(3),X(4))+Ux(t,X(1),X(2),X(3),X(4)),...
%          Ny(X(1),X(2), X(3),X(4))+Uy(t,X(1),X(2),X(3),X(4))]';
%        XP=@(t,X) [X(3),...
%            X(4),...
%            Linx(t),...
%            Liny(t)]';
%     XP =@(t,X) S1d(t);
    
     [T,S] = ode45(XP,[0,20*pi/b], S1S(0));
     plot(S(:,1),S(:,2),'b',L2(1),0,'g');
     hold;
     t1=linspace(0,2*pi/b,200);
     S2= zeros([4,200]);
     for k =1:200
         S2(:,k)=S1S(t1(k));
     end
     
     plot(S2(1,:),S2(2,:), 'r');

    
    
