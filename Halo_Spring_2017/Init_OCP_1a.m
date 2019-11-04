% PCR3B Time - Optimal Control Code Initialization Script

% Initial Conditions (Beginning State of Vehicle: [x, y, vx, vy]

global X0

x0 = -0.114720439757124;
y0 =  0.961641949132047;
vx0 = -0.0136629391404218;
vy0 =  0.118032720715187;

X0=[x0,y0,vx0,vy0]';

% Final Conditions (vehicle State at end of thrust)
global XF
% Position coordinates of (Halo or other) target orbit:
% xf = 0.461469276058880;
% yf = 0.866791033893503; 
% Vx = 0.000000377737727;
% Vy = -0.000000382078206;


xf =  0.722129624642116;
yf =  0.772097754592153;
Vx =  2.04719054208713e-016;
Vy = -0.135052580421232;

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
RBOUNDSQ = 1e6;



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
options1.Jacobian = 'finite';
% options1.Jacobian = 'broyden';
options1.Display = 'on';

% SOLUTION
% Lf =
% 
%          -1.30065764862467
%          0.022842957787406
%         -0.670761179318398
%         -0.772861038501534
%           1.75920606990477
% 
% 
% Fval =
% 
%     -2.33146835171283e-015
%     -7.32747196252603e-015
%     -8.03439346303003e-015
%      4.52415882534751e-015
%     -1.11022302462516e-014
% 
