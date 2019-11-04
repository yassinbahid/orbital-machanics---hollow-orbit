function Lf=OC_Run_Script(Initfile)
% Sample Commands to run the Optimal Control Software

% Set up and initializing data for a problem by
% Editing and running mfile script Init_OCP

disp(['Running ',Initfile,'  problem initialization script'])
eval(Initfile)

% Do a grid search to find good initial guesses for adjoints and final
% time: Lam0 = [lam1_0, lam2_0, lam3_0, lam4_0, tf]'
% Sample guess produced: Lam0 = [1 2 2 2 1]'
disp('Running  gridsearch_OC 5-D grid search to find good initial guesses for eqn solver mmfsolve')

[Lam0, dbest] = gridsearch_OC(.1,10,[-1,1],[-1,1],[-1,1],[-1,1],[2, 5.6727])

% Do a iterative search to refine an initial guesse for adjoints and final
% time: Lam0 = [lam1_0, lam2_0, lam3_0, lam4_0, tf]'
% Algorithm parameters for mmfsolve can be changed by setting options1
% structure.  See settings in Init_OCP
% Sample solution produced: Lf =[ 0.98176 1.9587 1.9593  1.9589   1 ]'
disp('Running  mmfsolve  iterative equation solver from initial guess')

[Lf,Fval,ExitFlag,Iters]=mmfsolve('shootingF_OC2',Lam0,options1)


% Plot Thrust Arc
PLOTON = 2;
d=shootingF_OC2(Lf);

% Plot Coast After Thrust Arc

global FINPT
tfc = 4;
Pf = Traj_simulation(tfc,FINPT(1:4));

% Plot Coast Without Thrust Arc
Pf = Traj_simulation(tfc,X0);

% Replot Trust Arc and Steering Angle
PLOTON = 2;
d=shootingF_OC2(Lf);
PLOTON = 0;