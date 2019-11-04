function [X0,VR VI] = init_halo(n,A,pair,theta)

% Function for producing state X0 in halo orbit about L_n
% (from solutions to linearized system Ydot = A Y, Y(0) = Y0
% Inputs:   Libration Point Number n
%           Initial conditions of diagonalized system: c = [c1  c2  c3  c4]
%           Notes on choosing c_k's  
%           1. Choose complex conjugate c_k's for corresponding complex
%           conjugate eigenvalues
%           2. Choose c_k = 0 for each eigenvalue whose mode
%           you wish to suppress.  
%           3. Choose |c| (modulus) smaller to bring initial
%           condition X0 closer to the Lagrange point.  This tends to improve
%           the accuracy of the resulting trajectory as an approx of the
%           nonlinear trajectory from the same X0.  
%           4. Choose complex phase angle (argument) of complex conjugate
%           c's to move the initial condition on its orbit.
%
% Outputs:  Numerical and Analytic Jacobian (J_n) Estimate at L_n
%           Eigenvalues and Eigenvectors of J_n
%           Initial Conditions X0 = L_n + Pc in nonlinear system
%           coordinates.  This can be used to start (or end) a PCR3BP

% Initialize Equilibria
% %Equilibria: (for mu=1/82.25)
% mu=1/82.25;
% Sample for finding Lagrange Points when mu is changed:
% x = newton_sysdalldamp('stateCR3B','stateCR3BJ',L(:,1),1.e-12,10,1,1)
% L(:,1) = [1.15571088526854      0 0 0]';
% L(:,2) = [0.83687838018159      0 0 0]'; 
% L(:,3) = [-1.00506575775433     0 0 0]';
% L(:,4) = [0.48784194528875   0.86602540378444   0 0]';
% L(:,5) = [0.48784194528875  -0.86602540378444   0 0]';

%Equilibria: (for mu=1/82.45)
global U
mu=U;
L(:,1) = [1.15559740258988      0 0 0]';
L(:,2) = [0.83702354452395      0 0 0]'; 
L(:,3) = [-1.00505347015937     0 0 0]';
L(:,4) = [0.48787143723469   0.86602540378444   0 0]';
L(:,5) = [0.48787143723469  -0.86602540378444   0 0]';

% Compute Jacobian of F(X) at Ln
% n=input('Enter number of Lagrange Point to Analyse.\n');

% Finite Difference Jacobian as a Check
% An = stateCR3BJ(L(:,n));
% An(1,3)=1;
% An(2,4) = 1;
% An(3,4)=2;
% An(4,3) = -2;
% An

% Analytic Jacobian 
% disp(['System Jacobian Matrix at L_',num2str(n)])
J = Jmatrix(L(1,n),L(2,n));

% Check on Analytic Jacobian 
% A-An
% A = An;  % Use F.F. Jacobian Until Analytic Partials are Fixed.



% Compute eigenvalues and eigenvectors
[P D]=eig(J);
% lam=zeros([4,1]);
% for k=1:4
% lam(k)=D(k,k);
% end
% disp('Eigenvalues of A')
% lam
% disp('Eigenvectors of A (columns)')
% P

% Set Initial Conditions
%c=input('Enter c_k''s to select term k of general solution.\n')';
if pair == 1|n == 1
    c = [A A 0 0]';
    b = imag(D(1,1));
    period = 2*pi/b;
    VR = real(P(:,1));
    VI = imag(P(:,1));
else
    c = [0 0 A A]';
    b = imag(D(3,3));
    period = 2*pi/b;
    VR = real(P(:,3));
    VI = imag(P(:,3));
end

% Convert c to U0
Y0=P*c;

%disp('Initial Conditions for Non-Linear System (in PCR3B Coordinates)')

t = (theta/(2*pi))*period;

X0 = L(:,n)+expm(J*t)*Y0;




