function [tf,ExitFlag] = ocp_time(X)
% returns optimal time in a 3-B control problem
% to the halo insertion point at angle alpha

if length(X) == 2
    alpha = X(2);
    r = X(1);
else
%     alpha = X;
%     r = .05;
    alpha = 2.3661;
    r = X;
end

% Set initial point for OCP
global X0
if isempty(X0)
X0 = circ_init(.3,pi,1);
end
% Compute final point for OCP
global XF
%r = .001;
XF = init_halo(4,r,1,alpha);

% Set initial guess for Newton iteration

global Lf

% Call Newton solver

options1 = mmfsolve('default');
options1.MaxIter = 70;
options1.Scale = 'on';
options1.FunTol = 1.e-9;
options1.Jacobian = 'finite';
%options1.Jacobian = 'broyden';
options1.Display = 'off';

[Lfn,Fval,ExitFlag,Iters]=mmfsolve('shootingF_OC2',Lf,options1);

% Check if Solved
if ExitFlag<=1
    Lf = Lfn;
end
% Unload final time

tf = Lf(5);