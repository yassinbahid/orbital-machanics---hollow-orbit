% Program by Bruce Lundberg, PHD, and Yassin Bahid
% PCR3B Time - Optimal Control Code Initialization Script

% Initial Conditions (Beginning State of Vehicle: [x, y, vx, vy]
global X0
x0 = S10(1);
y0 = S10(2);
vx0 = S10(3);
vy0 = S10(4);

X0=[x0,y0,vx0,vy0]';

% Final Conditions (vehicle State at end of thrust)
global XF
% Position coordinates of (Halo or other) target orbit:
% xf = 0.461469276058880;
% yf = 0.866791033893503; 
% Vx = 0.000000377737727;
% Vy = -0.000000382078206;
% xf =  0.71927; 
% yf =  0.48576;
% Vx = -0.49203;
% Vy = -0.41552;

xf =  x0;
yf =  y0;
Vx = vx0;
Vy = vy0;

XF = [xf,yf,Vx,Vy]';

% Global Thrust Variable and Value

% Reduced Mass Parameter mu = m1/(m1+m2)
global U
U = mu;




%____________________________________________________________________

% ALGORITHM PARAMETERS
%____________________________________________________________________

% Trajectory Plot Switch (= 0 for off)
global PLOTON
PLOTON = 0;

% Set Integrator Options

% Numerical Integrator Relative and Absolute Error Tolerances
global TOLR TOLA
TOLR = 1.e-12;
TOLA = 1.e-12;

% Set Escape Event Parameter RBOUNDSQ (to stop integration if ||X(t)|| gets
% too big
global RBOUNDSQ
RBOUNDSQ = 1e6;

% Set Iterative Nonlinear Equation Solver Options
options1 = mmfsolve('default');
options1.MaxIter = 50;
options1.Scale = 'on';
options1.FunTol = 1.e-8;
options1.Jacobian = 'finite';
% options1.Jacobian = 'broyden';
options1.Display = 'on';