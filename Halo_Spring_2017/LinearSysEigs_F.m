function [X,V_R,V_I,b,L]=Linear_SysEigs_F(r,alpha,n_lib,n_eig)

% Script for analysing solutions to linear system Ydot = A Y, Y(0) = Y0
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

switch n_lib
    case 1
        L = [1.15559740258988      0 0 0]';
    case 2
        L = [0.83702354452395      0 0 0]';
    case 3
        L = [-1.00505347015937     0 0 0]';
    case 4
        L = [0.48787143723469   0.86602540378444   0 0]';
    case 5
        L = [0.48787143723469  -0.86602540378444   0 0]';
    otherwise
        disp('n_lib must be a libration point number, 1, ..., 5')
end

% Compute Jacobian of F(X) at Ln
% Finite Difference Jacobian as a Check
% An = stateCR3BJ(L(:,n));
% An(1,3)=1;
% An(2,4) = 1;
% An(3,4)=2;
% An(4,3) = -2;
% An

% % Analytic Jacobian 
% disp(['System Jacobian Matrix at L_',num2str(n)])
 A = Jmatrix(L(1),L(2));

% Compute eigenvalues and eigenvectors
[P D]=eig(A);
lam=zeros([4,1]);
for k=1:4
lam(k)=D(k,k);
end
% disp('Eigenvalues of A')
% lam
% disp('Eigenvectors of A (columns)')
% P
% 
% % Set Initial Conditions
% c=input('Enter c_k''s to select term k of general solution.\n')';
cc = r*exp(i*alpha);
switch n_eig
    case 1
        c = [cc, conj(cc), 0, 0]';
        V_R = real(P(:,1));
        V_I = imag(P(:,1));
        b = imag(D(1,1));
    case 2
        c = [0, 0, cc, conj(cc)]';
        V_R = real(P(:,3));
        V_I = imag(P(:,3));
        b = imag(D(3,3));
    otherwise
        disp('n_eig must be an eigenvalue pair number 1 or 2')
end

% Convert c to U0
Y=P*c;

% disp('Initial Conditions for Non-Linear System (in PCR3B Coordinates)')
X = L+Y;


% lammin = realmax;
% for k=1:4   
%     if ~isreal(lam(k))
%         lamm = abs(imag(lam(k)));
%         if lamm<lammin
%             lammin=lamm;
%         end
%     end
% end
% period = 2*pi/lammin;
% 
% % Compute Trajectory of Linear System
% 
% t = 0:.001*period:2*period;
% y = zeros([4,length(t)]);
% for k = 1:length(t)
%     y(:,k) = expm(A*t(k))*Y0 + L(:,n);
% end
% 
% 
% 
% 
% subplot(2,1,1)
% plot(y(1,:),y(2,:),'b',y(1,1),y(2,1),'rp')
% title(['Position:  PCR3B Linearized about L_',num2str(n)])
% xlabel('x-axis')
% ylabel('y-axis')
% hold on
% 
% axis equal
% % v=axis;
% % axis(1.1*v)
% 
% 
% subplot(2,1,2)
% plot(y(3,:),y(4,:),'b',y(3,1),y(4,1),'rp')
% title(['Velocity:  PCR3B Linearized about L_',num2str(n)])
% xlabel('V_x-axis')
% ylabel('V_y-axis')
% hold on
% 
% 
% axis equal
% v=axis;
% axis(1.1*v)
% 
% % Coordinate axes (eigenvectors)
% iflag = 0;
% for k =1:4
%     
%     % skipping plot?   
%     if imag(lam(k))> 0
%         continue
%     end
%     
%     skip=input(['To skip plot of eigenvector ',num2str(k),' enter s  '],'s');
%     if skip=='s'
%         continue
%     end
%     
%     YR=.5*norm(c)*real(P(1:2,k))/norm(real(P(1:2,k)));
%     VR=.5*norm(c)*real(P(3:4,k))/norm(real(P(3:4,k)));
%     subplot(2,1,2)
%     quiver(0,0,VR(1),VR(2))
%     subplot(2,1,1)
%     quiver(L(1,n),L(2,n),YR(1),YR(2))
%     
%     if any(imag(P(:,k)))
%         iflag = 1;
%         YI=.5*norm(c)*imag(P(1:2,k))/norm(imag(P(1:2,k)));
%         VI=.5*norm(c)*imag(P(3:4,k))/norm(imag(P(3:4,k)));
%         quiver(L(1,n),L(2,n),YI(1),YI(2))
%         subplot(2,1,2)
%         quiver(0,0,VI(1),VI(2))
%     end 
%        
% end

% % Function to Compute the Analytic Jacobian
% function A=Amatrix(x,y,u)
% % A Matrix: Project, Lundberg, James, Cole
% 
% 
% % General Parts
% mus=(1-u);
% r1s=((x+u).^2+y.^2);
% r2s=((x-1+u).^2+y.^2);
% 
% % Some Partials
% F_3x = 1 - mus*(y.^2 -2*(x + u).^2)./(r1s.^2.5)  -  u*(y.^2 -2*(x + u-1).^2)./(r2s.^2.5);
% 
% F_3y =    3*mus*((x + u).*y)./(r1s.^2.5)  +  3*u*((x + u-1).*y)./(r2s.^2.5);
% 
% F_4x =    3*mus*((x + u).*y)./(r1s.^2.5)  +  3*u*((x + u-1).*y)./(r2s.^2.5);
% 
% F_4y = 1 - mus*((x + u).^2 - 2*y.^2 )./(r1s.^2.5)  -  u*((x + u-1).^2-2*y.^2)./(r2s.^2.5);
% 
% 
% % The Jacobian of F for the PCR3BP Xdot = F(X)
% A=[0,0,1,0;0,0,0,1;F_3x,F_3y,0,2;F_4x,F_4y,-2,0];
