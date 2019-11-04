% PCR3B Time - Optimal Control Code Initialization Script

% Initial Conditions (Beginning State of Vehicle: [x, y, vx, vy]

global X0

% x0 = -0.114720439757124;
% y0 =  0.961641949132047;
% vx0 = -0.0136629391404218;
% vy0 =  0.118032720715187;

x0 =        -0.262128562765313; % r = .4, th = 4pi/3, dir = -1
y0 =        -0.433012701892219;
vx0 =          -1.6503077178403;
vy0 =         0.952805605140814;


X0=[x0,y0,vx0,vy0]';

% Final Conditions (vehicle State at end of thrust)
global XF

% Position coordinates of (Halo or other) target orbit:

xf =        -0.412128562765313;  %r = .8, th = 4pi/3, dir= -1
yf =        -0.692820323027551;
Vx =       -1.65517653171928;
Vy =        0.955616616144479;
% 
% xf =         0.237871437234688; %r = .5, th = pi/3, dir = -1
% yf =         0.433012701892219;
% Vx =           1.6503077178403;
% Vy =        -0.952805605140813;
% 
% xf =        -0.262128562765312; %r = .5, th = 2pi/3, dir = -1
% yf =         0.433012701892219;
% Vx =           1.6503077178403;
% Vy =         0.952805605140812;


XF = [xf,yf,Vx,Vy]';

% Global Thrust Variable and Value

global TH
TH = 1;


% Reduced Mass Parameter mu = m1/(m1+m2)
global U
U = 1/82.45;

% Set Escape Event Parameter RBOUNDSQ (to stop integration if ||X(t)|| gets
% too big
global RBOUNDSQ
RBOUNDSQ = 1e8;

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


% Set Iterative Nonlinear Equation Solver Options

options1 = mmfsolve('default');
options1.MaxIter = 50;
options1.Scale = 'on';
options1.FunTol = 1.e-12;
%options1.Jacobian = 'finite';
 options1.Jacobian = 'broyden';
options1.Display = 'on';

% Results
% to r = .5, th = pi/3
% Lf =
%         -0.238809484076768
%           2.37185280355471
%        -0.0542063383321254
%          0.291092348820945
%          0.812474557121304
 
