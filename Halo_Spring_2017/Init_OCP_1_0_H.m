% PCR3B Time - Optimal Control Code Initialization Script

% Initial Conditions (Beginning State of Vehicle: [x, y, vx, vy]

global X0

% x0 = -0.114720439757124;
% y0 =  0.961641949132047;
% vx0 = -0.0136629391404218;
% vy0 =  0.118032720715187;

% x0 =        -0.262128562765313;
% y0 =        -0.433012701892219;
% vx0 =          -1.6503077178403;
% vy0 =         0.952805605140814;
% 
% x0 =        -0.412128562765313;  %r = .8, th = 4pi/3, dir= -1
% y0 =        -0.692820323027551;
% vx0 =       -1.65517653171928;
% vy0 =        0.955616616144479;

% x0 =        -0.312128562765312; %r = .3, th = pi, dir= 1
% y0 =     0;
% vx0 =    0;
% vy0 =          -1.5146362695544;
% 
% X0=[x0,y0,vx0,vy0]';

X0 =[-0.168128562765312; 0;0 ; -2.3604477369981];


% Final Conditions (vehicle State at end of thrust)
global XF
% Position coordinates of (Halo or other) target orbit:
% xf = 0.461469276058880;
% yf = 0.866791033893503; 
% Vx = 0.000000377737727;
% Vy = -0.000000382078206;

xf = 0.476158527864319; 
yf = 0.863647957647671;
Vx = -0.0111810642084073;
Vy =  0.0112357640106497;

XF = [xf,yf,Vx,Vy]';
[XF,VR VI] = init_halo(4,.05,1,pi/1.67);

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
RBOUNDSQ = 1e-1;
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
options1.FunTol = 1.e-9;
options1.Jacobian = 'finite';
%options1.Jacobian = 'broyden';
options1.Display = 'on';

% SOLUTIONS 
% x0 =        -0.312128562765312; %r = .3, th = pi, dir= 1  ... 
% X0 = circ_init(.3,pi);
% Lf =[         -2.44721447606011
%          -1.07377916385952
%         -0.272395673952003
%         -0.471701260583287
%           1.69022110889836];
% x0 =        %r = .21, th = pi, dir= 1  ... 
% Lf =[  155.9698
%    -6.8083
%    -0.0333
%    13.7282
%     2.1828]

% x0 =        %r = .156, th = pi, dir= 1  ... 
%  Solution by parameter tracking from r = .3 solution
% Lf =[101.349
%          -3.17729
%         -0.0704
%         5.9687
%           5.6727];
% Solution using Xbest from  
% [Xbest, dbest] = gridsearch_OC3([1,1],0,[-1,1],15,[-1,1],15,[-1,1],15,0:.02:5,1)
% Lfn = [0.63860440712209
%         -0.935088268049324
%        -0.0861929529572342
%        -0.0959953421263096
%           2.24494590758726];
% x0 =        %r = .1, th = pi, dir= 1  ...
%  [Xbest, dbest] = gridsearch_OC3([1,1],0,[-1,1],19,[-1,1],19,[-1,1],19,0:.01:4,1)
% X0=circ_init(.1,pi);
% Lfn =[ -11.2680316231566
%           1.61047583776428
%         0.0463155725590549
%         -0.385253561001683
%           3.63132125482941];

% [XF,VR VI] = init_halo(4,.05,1,pi/1.67); and 
% x0 =        -0.312128562765312; %r = .3, th = pi, dir= 1  ... 
% Lfn =[         -2.48143916083695
%          -1.11297616234209
%         -0.280290088764464
%         -0.465932454931842
%           1.66968519273681]
