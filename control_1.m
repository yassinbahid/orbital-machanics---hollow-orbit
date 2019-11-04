function Uc = control_1(t,X)
% function Uc = control_1(t,X)
% Computes the "feedback control" needed to compensate for differences
% between linear and nonlinear accelerations.
 global J Y0 
% r = R;
% n_lib = NLIB;
% n_eig=NEIG;
% 
% 
% [X0,V_R,V_I,b,A,Y0]=LinearSysEigs_F(r,0,n_lib,n_eig);
% period = 2*pi/b;
% alpha = t*b;
% [X,V_R,V_I,b,A,Y]=LinearSysEigs_F(r,alpha,n_lib,n_eig);
[X1,SL] =Lin_Traj(t,Y0);
S_dot=state_dot(t,X);
%S_dot=state_dot(t,X);
Uc = J*SL - S_dot;