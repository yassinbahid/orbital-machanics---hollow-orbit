function XP = PCR3BC(t, X)
% Evaluates RHS of the equations of motion
% for the planar circular restricted three-body problem
% with Lyapunov feedback control using final state Xf = (xf, yf, uf, vf)
% and global gain matrix (2 x 4) GAINS
% and trust/mass = TH, and primaries mass ratio U
% as initialized in the calling workspace
% TH_FLAG = 
%           -1  Steepest Descent of J = Jacobi Constant of the Motion
%            1  Steepest Ascent of J = Jacobi Constant of the Motion
%           -2  Steepest Descent of Lyapunov Function J = 0.5(X - XF)^T*G*(X - XF)
%            2  Steepest Ascent of Lyapunov Function  J = 0.5(X - XF)^T*G*(X - XF)


% Fixed global parameters
global TH U TH_FLAG
global GAINS XF

% Unload state components
x = X(1);
y = X(2);
u = X(3);
v = X(4);

% Steepest Descent (or Ascent) Controls from Lyapunov Descent on J = Jacobi
% Constant or V = (X(t)-XF)^T*GAINS*(X(t)-XF)
if TH_FLAG == -1
    T1 = -u;
    T2 = -v;
elseif TH_FLAG == 1
    T1 = u;
    T2 = v;
elseif TH_FLAG == -2
    CNTRL = -GAINS*(X - XF);
    T1 = CNTRL(1);
    T2 = CNTRL(2);
elseif TH_FLAG == 2
    CNTRL = GAINS*(X - XF);
    T1 = CNTRL(1);
    T2 = CNTRL(2);
end

sqrt12 = sqrt(T1.*T1+T2.*T2);
if sqrt12 > 0
   T1 = TH*T1./sqrt12;
   T2 = TH*T2./sqrt12;
end

% State Equations of motion
r13 =  sqrt((x-( -U)).^2 + y.^2).^3;
r23 =  sqrt((x-(1-U)).^2 + y.^2).^3;
XP(1) = u;
XP(2) = v;
XP(3) = x + 2*v - (1-U)*(x+U)./r13 - U*(x-(1-U))./r23 + T1;
XP(4) = y - 2*u - (1-U)*y./r13 - U*y./r23 + T2;

XP=XP';





