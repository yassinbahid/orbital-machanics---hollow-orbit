function S_dot=state_dot_OC2(t,S)
% Program by Bruce Lundberg, PHD, and Yassin Bahid
% This function will take the 4-input compnents described in the State(S)
% and find the rates of change for each of those components in the form of
% ODE's. The specific form of the ODE's adhers to the
% circular-plane-restricted 3-body problem.

% In addition the function will construct 4-additional ODE's (adjoints)
% which will be the four components of the Adjoint equation constraining
% our optimal control problem. Adjoint equation is:
%                               
%                           Lambda_dot = - Hamiltonian(partial with
%                                        respect to state S)Transpose 
                                                                                     
% Necessary paramters:
global U
u=U;

% Load state variables into locals
      x = S(1);
      y = S(2);
      vx = S(3);
      vy = S(4);

% Load adjoint variables into locals
      lam1= S(5);
      lam2= S(6);
      lam3= S(7);
      lam4= S(8);
      

% Compute local quantities
      um1 = 1-u;

      xpu = x + u;
      xpum1 = x + u - 1;
      xpus = xpu*xpu;
      xpum1s = xpum1*xpum1;
      ys = y*y;

      r1s = xpus + ys;
      r2s = xpum1s + ys;
      r1 = sqrt(r1s);
      r2 = sqrt(r2s);

      r13=r1s*r1;
      r23=r2s*r2;
      r15=r1s*r13;
      r25=r2s*r23;


% Compute thrust angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE CHANGE CHANGE CHANGE CHANGE CHANGE CHANGE CHANGE CHANGE CHANGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     pheta=datan2(-lam4,-lam3)
%       l34n = sqrt(lam3*lam3+lam4*lam4);
%       cth = -lam3/l34n;
%       sth = -lam4/l34n;


% Formulation of State ODEs for x, y xdot, ydot:
      um1or13 = um1/r13;
      uor23 = u/r23;
      S_dot(1) = vx;
      S_dot(2) = vy;
      S_dot(3) = 2*vy + x - um1or13*xpu - uor23*xpum1-.5*lam3;
      S_dot(4) =-2*vx + y - y*(um1or13   + uor23)     -.5*lam4;

% Adjoint ODEs:   dlam = -DF^T lam
      
      um1or15 = um1/r15;
      uor25 = u/r25;
      F31 = 1 - um1or15*(r1s - 3*xpus) - uor25*(r2s - 3*xpum1s);
      F32 = 3*y*(um1or15*xpu + uor25*xpum1);
      F41 = F32;
      F42 = 1 - um1or15*(r1s - 3*ys) - uor25*(r2s - 3*ys);

      S_dot(5) =                - F31*lam3 - F41*lam4;
      S_dot(6) =                - F32*lam3 - F42*lam4;
      S_dot(7) = -lam1                     +   2*lam4;
      S_dot(8) =        -lam2   - 2*lam3;


S_dot=S_dot';   % stored as a column vector.


