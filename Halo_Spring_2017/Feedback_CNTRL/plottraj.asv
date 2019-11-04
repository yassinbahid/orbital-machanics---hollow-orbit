function plottraj(t1, t2,X0,Xf,Gains,thrustI,thrustC,thrust_ang,orbline,thrstwgt2)
% integrates Planar circular restricted 3-Body ode's 
% with Lyapunov control with thrust magnitude Thrust
% and gains Gains=[G1, ... , G4];
% if Gains = [  ] then control based on minimizing Jacobi Integral used
% Sample Commands:  X0=[1.2; 0; 0; -1.04935750983031990726];
%          plottraj(10,0,X0,Xf,Gains,0,0,0,'m-')   % Coast in periodic orbit
%          plottraj(30,0,X0,Xf,Gains,0,.01,0,'m-') % Flys near L2 and escapes
% XfL2 = [L(2),0,0,0]';,plottraj(12.4152884,10,X0,XfL2,0,.9871,.003793249,pi/2-.18500,'m-',0)
%          plottraj(12.415,0,X0,XfL2,0,.987,.00378,pi/2,'m-')  Flys to L2
%          plottraj(12.415,20,X0,XfL2,0,.987,.00378,pi/2,'m-',0) and turn
%          off low thrust
%  XfL2 = [L(2),0,0,0]';,plottraj(12.4152884,10,X0,XfL2,0,.98703,.003793249,pi/2-.184710,'m-',0)
%  >> GAINS=[1 0 7 1;0 1 1 7],    plottraj(1.91,4,X0,XF,[ ] ,0,1,0,'m-',0)
%     X0=[1.2; 0; 0; -1.04935750983031990726]';
%     Xf = [L(2),0,0,0]';,
%
% >>Init_OCP_1_0_H
% >> plottraj(.451,12.4,X0,XF,[ ],0,1,0,'g-',0)
% >>% >>Gains = [5 1 50 1;1 5 1 20];      
% >>plottraj(.647,15.7,X0,XF,Gains,0,1,0,'c-',0)
   

% Declare global parameters

global xf yf uf vf
global TH U TH_FLAG TH_SAVE
global GAINS XF L

TH_FLAG = 1;   % STEEPEST DESCENT OF JACOBI POTENTIAL  (+1 for ASCENT)
if ~isempty(Gains)
    GAINS = Gains;
    TH_FLAG = -2;  % CHANGE TO 2 FOR ESCAPE FROM XF
    TH_FLAG = 2;  % CHANGE TO -2 FOR ESCAPE FROM XF

    XF = Xf;
end
    

% Initialize Primaries Mass Ratio and target final state
% MU = 1/82.45;   % Earth-Moon System
xf = Xf(1);
yf = Xf(2);
uf = Xf(3);
vf = Xf(4);

TH_SAVE = TH;
TH = thrustC;
TI =thrustI;

% Set ode45 algorthm parameters
options1 = odeset('AbsTol', 1.e-12,'RelTol', 1.e-12);

if t1>0
% Phase 1
Tspan1 = [0, t1];
[tout1,Xout1]=ode45('PCR3BC',Tspan1,X0,options1);
else
    Xout1 = X0;
    tout1 = 0;
end

% Apply Impulse
ctheta = cos(thrust_ang);
stheta = sin(thrust_ang);
X02=Xout1(end,:)+[0,0,TI*ctheta,TI*stheta];

if t2>0        
%  Phase 2 
TH = thrstwgt2*TH;
Tspan2 = [t1, t1+t2];
[tout2,Xout2]=ode45('PCR3BC',Tspan2,X02,options1);


% Combine and Unload Components

Xout =[Xout1;Xout2(1:end,:)];
tout = [tout1;tout2(1:end)];

else
   Xout =[Xout1;X02];
   tout = [tout1;tout1(end)];
end

x = Xout(:,1);
y = Xout(:,2);
u = Xout(:,3);
v = Xout(:,4);

% Plotting Bodies and Libration Points

L=eulerpts(U);  % computes Euler colinear points

%subplot(2,1,1)
plot(-U,0,'bo')
hold on
plot(1-U,0,'co')
L(4) = .5-U;
L(5) = sqrt(3)/2;

plot(L(1), 0,'b.')
plot(L(2), 0,'m.')
plot(L(3), 0,'c.')
plot(L(4),-L(5),'r.')
plot(L(4), L(5),'r.')
plot(xf,yf,'gs')

% Plot Trajectory
plot(x(1),y(1),'rp')
plot(x,y,orbline)
legend('Earth','Moon','L_1','L_2','L_3', 'L_4',' L_5', 'Target','X_0','X(t)',3)
title('Relative Trajectory')
xlabel('x')
ylabel('y')
% Plot Velocities

figure
%subplot(2,1,2)
plot(u(1),v(1),'rp')
hold on
plot(u,v,orbline)
title('Relative Velocity')
xlabel('v_x')
ylabel('v_y')
grid on

figure
%subplot(2,1,2)
r1 =  sqrt((x-( -U)).^2 + y.^2);
r2 =  sqrt((x-(1-U)).^2 + y.^2);
Jout= 0.5*(u.^2+v.^2 - x.^2 - y.^2);
Jout = Jout-(1-U)./r1 - U./r2;

% Compute Jacobi Constant for L1
rL11 =  sqrt((L(1)-( -U)).^2 + 0.^2);
rL12 =  sqrt((L(1)-(1-U)).^2 + 0.^2);
JoutL1= 0.5*(0.^2+0.^2 - L(1).^2 - 0.^2);
JoutL1 = JoutL1-(1-U)./rL11 - U./rL12;
% Compute Jacobi Constant for L2
rL21 =  sqrt((L(2)-( -U)).^2 + 0.^2);
rL22 =  sqrt((L(2)-(1-U)).^2 + 0.^2);
JoutL2= 0.5*(0.^2+0.^2 - L(2).^2 - 0.^2);
JoutL2 = JoutL2-(1-U)./rL21 - U./rL22;
% Compute Jacobi Constant for L3
rL31 =  sqrt((L(3)-( -U)).^2 + 0.^2);
rL32 =  sqrt((L(3)-(1-U)).^2 + 0.^2);
JoutL3= 0.5*(0.^2+0.^2 - L(3).^2 - 0.^2);
JoutL3 = JoutL3-(1-U)./rL31 - U./rL32;

% Compute Jacobi Constant for L4 and L5
rL51 =  sqrt((L(4)-( -U)).^2 + L(5).^2);
rL52 =  sqrt((L(4)-(1-U)).^2 + L(5).^2);
JoutL5= 0.5*(0.^2+0.^2 - L(4).^2 - L(5).^2);
JoutL5 = JoutL5-(1-U)./rL51 - U./rL52;
%plot(tout,Jout-Jout(1),orbline)
tJ=tout(1);
plot(tout,Jout,orbline,tJ,JoutL1,'b.',tJ,JoutL2,'m.',tJ,JoutL3,'c*',tJ,JoutL5,'r*')
title('Jacobi Integral (Conserved Quantity)')
xlabel('t')
ylabel('J')
% legend('J(t)','J(L_1)','J(L_2)','J(L_3)', 'J(L_4 and L_5)','Location','SouthWest')
legend('J(t)','J(L_1)','J(L_2)','J(L_3)', 'J(L_4 and L_5)',3)
grid on
axis([0,tout(end),-4,0])

% Return Thrust to Original Value
TH = TH_SAVE;