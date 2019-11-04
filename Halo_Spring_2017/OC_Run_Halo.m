function Lf=OC_Run_Halo(Initfile)
% Sample Commands to run the Optimal Control Software

% Set up and initializing data for a problem by
% Editing and running mfile script Init_OCP

disp(['Running ',Initfile,'  problem initialization script'])
eval(Initfile)

% Do a grid search to find good initial guesses for adjoints and final
% time: Lam0 = [lam1_0, lam2_0, lam3_0, lam4_0, tf]'
% Sample guess produced: Lam0 = [1 2 2 2 1]'
disp('Running  gridsearch_OC 5-D grid search to find good initial guesses for eqn solver mmfsolve')
global TF
%[Lam0, dbest] = gridsearch_OC_Halo([-1,1],4,[-1,1],4,[-1,1],4,[-1,1],4,[TF,TF+.1],1)
[Lam0, dbest] = gridsearch_OC_Halo([1,1],1,[0,0],1,[0,0],1,[1,1],1,[TF,TF],1)

% Do a iterative search to refine an initial guesse for adjoints and final
% time: Lam0 = [lam1_0, lam2_0, lam3_0, lam4_0, tf]'
% Algorithm parameters for mmfsolve can be changed by setting options1
% structure.  See settings in Init_OCP
% Sample solution produced: Lf =[ 0.98176 1.9587 1.9593  1.9589   1 ]'
disp('Running  mmfsolve  iterative equation solver from initial guess')

%Lam0=[0.0007   -0.0060   -0.0039   -0.0060]'; % Initial guess lam = -2*U(0,X0) Feedback
Lam0=1.0e-004*[   -0.4753     0.0751   -0.1240    -0.0570]'; % OCP3 Solution for r=.001
% Lam0=1.0e-003*[    0.3013   -0.3229    0.1446   -0.1915]'; % OCP4 Solution for r=.001

[Lf,Fval,ExitFlag,Iters]=mmfsolve('shootingF_OC2H',Lam0,options1)


% Plot Thrust Arc
PLOTON = 2;
d=shootingF_OC2H(Lf);

% Plot Coast After Thrust Arc

global FINPT
tfc = 4;
Pf = Traj_simulation(tfc,FINPT(1:4));

% Plot Coast Without Thrust Arc
Pf = Traj_simulation(tfc,X0);

% Replot Trust Arc and Steering Angle
PLOTON = 2;
d=shootingF_OC2H(Lf);
PLOTON = 0;