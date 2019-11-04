function S_dot=state_dot_Equ(S)

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

u1 = 1-u; 
r1 = sqrt((S(1)+u)^2+S(2)^2);
r2 = sqrt((S(1)-u1)^2+S(2)^2);
% lam1= S(5); % come in as the last four of eight components in S.
% lam2= S(6);
% lam3= S(7);
% lam4= S(8); 
%    

% pheta=atan2(-lam4,-lam3); % Thrust angle.
% global TH;
% TH_local=TH;

% Formulation of ODE's:

S_dot(1) = S(3);
S_dot(2) = S(4);                                          
S_dot(3) = 2*S(4)+S(1)-((u1*(S(1)+u))/(r1^3))-((u*(S(1)-u1))/(r2^3)); %+ TH_local*cos(pheta); % Equals second derivative of x-position 
S_dot(4) = -2*S(3)+S(2)-((u1*S(2))/(r1^3))-((u*S(2))/(r2^3)); %+ TH_local*sin(pheta); % Equals second derivative of y-position 

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


