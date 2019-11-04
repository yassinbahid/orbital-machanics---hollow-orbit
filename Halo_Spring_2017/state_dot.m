function S_dot=state_dot(t,S)

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

global U;
u=U;

% Formulation of ODE's:
% Load state variables into locals
      x = S(1);
      y = S(2);
      vx = S(3);
      vy = S(4);

% % Load adjoint variables into locals
%       lam1= S(5);
%       lam2= S(6);
%       lam3= S(7);
%       lam4= S(8);

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
%       r15=r1s*r13;
%       r25=r2s*r23;


% % Compute thrust angle
% %     pheta=datan2(-lam4,-lam3)
%       l34n = sqrt(lam3*lam3+lam4*lam4);
%       cth = -lam3/l34n;
%       sth = -lam4/l34n;


% Formulation of State ODEs for x, y xdot, ydot:
      um1or13 = um1/r13;
      uor23 = u/r23;
      S_dot(1) = vx;
      S_dot(2) = vy;
      S_dot(3) = 2*vy + x - um1or13*xpu - uor23*xpum1; %+ Thrust*cth;
      S_dot(4) =-2*vx + y - y*(um1or13   + uor23);     %+ Thrust*sth;


% S_dot(1) = S(3);
% S_dot(2) = S(4);                                          
% S_dot(3) = 2*S(4)+S(1)-((u1*(S(1)+u))/(r1^3))-((u*(S(1)-u1))/(r2^3)); %+ TH_local*cos(pheta); % Equals second derivative of x-position 
% S_dot(4) = -2*S(3)+S(2)-((u1*S(2))/(r1^3))-((u*S(2))/(r2^3)); %+ TH_local*sin(pheta); % Equals second derivative of y-position 

% % Adjoints:
% 
% S_dot(5) = ( 1-((((S(1)-(1-u))^2+S(2)^2)*u-3*(u*S(1)+u)*(S(1)-(1-u)))/(r1^5)) - (((1-u)*((S(1)+u)^2+S(2)^2)^(1/2)*(-2*(S(1)+u)^2+S(2)^2))/(r1^6)) )*lam3 + ( -((-3*S(2)*(1-u)*(S(1)+u))/(r2^5)) - ((-3*u*S(2)*(S(1)-1+u))/(r1^5)) )*lam4;
% S_dot(5) = -S_dot(5);
% S_dot(6) = ( (3*S(2)*(u*S(1)+u))/(r1^5) + (3*S(2)*(1-u)*(S(1)+u))/(r2^5) )*lam3 + ( -((1-u)*((r2^3)-3*S(2)^2*r1))/(r1^6) - ((S(1)-(1-u))^2-2*S(2)^2)/(r1^5) )*lam4;
% S_dot(6) = -S_dot(6);
% S_dot(7) = lam1 - 2*lam4;
% S_dot(7) = -S_dot(7);
% S_dot(8) = lam2 + 2*lam3;
% S_dot(8) = -S_dot(8);



S_dot=S_dot';   % stored as a column vector.


