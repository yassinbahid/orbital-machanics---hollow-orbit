function [d,TE] = shootingF_OC2(X)

% Program shootingF_new.m. This program will integrate the system of 8 differential
% equations defined in the function statedot_v2.m forward in time according
% to a shooting method where a set of 5 initial parameters are described in X.  When the ouput 
% error d approaches the zero vector, the parameters described in X approach 
% a viable solution.

% X=initial lambda values plus final time.
% X0=initial state of the ship.

% ATTENTION !!! In order to accomodate the thrust term in the ode system
%               setup the user must declare; global TH, and enter a value for TH.
% (You can run the script Init_OCP to do this)



% Sample Input:

% Trial 1:
% global TH , TH=0;
% X0=[0.244647414873717,0.836089636635792,1.142126232466290,0.489143510530548];
% d = shootingF_adj([0,.1,.1,-1,.1],X0)

% Trial 2: 
% global TH  , TH=0;
% X0=[0.244647414873717,0.836089636635792,1.142126232466290,0.489143510530548];
% d = shootingF_adj([-0.2000   -0.2000   -0.2000   -0.2000    1.0000],X0)

% Unloading intial parameters:
lam1=X(1);
lam2=X(2);
lam3=X(3);
lam4=X(4);
tf=X(5);

%Lam1=lam1;
%Lam2=lam2;
%Lam3=lam3;
%Lam4=lam4;
%TF=tf;
%global Lam1 Lam2 Lam3 Lam4 TF

% Initial state of the ship:
%X0=[0.244647414873717,0.836089636635792,1.142126232466290,0.489143510530548];
global X0
% Load Initial Conditions for State Adjoint Integration
M0 = [X0;X(1:4)];

% Load Final State of ship components
global XF
xf = XF(1);
yf = XF(2); 
Vx = XF(3);
Vy = XF(4);


% Integration Forward:
global TOLR TOLA
global TF
global U

options1 = odeset('RelTol',TOLR,'AbsTol',TOLA,'Events', @events);
[T,M,TE,IE]=ode45(@state_dot_OC,[0,tf],M0,options1);
Mf = M(end,:)';

if ~isempty(TE)
    disp(['Func Eval ERROR:  TE = ',num2str(TE)])
end

% final state; after integration:
x_tf=Mf(1);
y_tf=Mf(2);
u_tf=Mf(3);  % x-component of velocity.
v_tf=Mf(4);  % y-component of velocity.
lam1_tf=Mf(5);
lam2_tf=Mf(6);
lam3_tf=Mf(7);
lam4_tf=Mf(8);

% Calculation of Error from desired target Halo. 
d(1) = x_tf-xf;
d(2) = y_tf-yf;
d(3) = u_tf-Vx;
d(4) = v_tf-Vy;

% Hamiltonian/Cap_phi Transversality condition: 
f_tf = state_dot_OC(tf,M(end,:)');    % need u_dot and v_dot at tf.
f_tf = f_tf(1:4);

%d(5) = 1+ lam1_tf*u_tf + lam2_tf*v_tf + lam3_tf*f_tf(3) + lam4_tf*f_tf(4);
d(5) = 1+ [lam1_tf , lam2_tf , lam3_tf , lam4_tf]*f_tf;



% Plotting Section:
global PLOTON
if PLOTON ~= 0     % In case ~= 0, plots will be produced for reference and refinement.
    subplot(2,1,1)
    plot(M(1,1),M(1,2),'rp')  % initial position.
    hold on
    plot(Mf(1),Mf(2),'gp')    % final point
    plot(xf,yf,'bo')  % targeted position
    plot(M(:,1),M(:,2))       % trajectory.
%    legend('X_0','X(t_f) Attained', 'X_f Target')
    xlabel('x-axis'), ylabel('y-axis'), title('xy-trajectory plot')
    
    subplot(2,1,2)
    if PLOTON == 1
%         if x_tf-X0(1) > 1  &   y_tf-X0(2) > 1
%             text(x_tf-.1*x_tf, y_tf-.1*y_tf,['(',num2str(x_tf),' ',num2str(y_tf),') ',num2str(d(3)), '  ', num2str(d(4)),'  ' num2str(d(5))])
%             axis equal
%         else
%             avg_x = abs(mean( x_tf-X0(1)));
%             avg_y = abs(mean( y_tf-X0(2)));
%             text(avg_x,avg_y,['(',num2str(x_tf),' ',num2str(y_tf),') ',num2str(d(3)), '  ', num2str(d(4)),'  ' num2str(d(5))])
%         end
       
        plot(M(1,3),M(1,4),'rp')  % initial velocity.
        hold on
        plot(Mf(3),Mf(4),'gp')    % final velocity
        plot(Vx,Vy,'bo')  % targeted velocity
        plot(M(:,3),M(:,4))       % trajectory.
%        legend('V_0','V(t_f) Attained', 'V_f Target')
        xlabel('Vx-axis'), ylabel('Vy-axis'), title('Vx Vy-trajectory plot')
        
    else

        pheta=atan2(-M(:,8),-M(:,7)); % Thrust angle.
        cut = -2;
        ibn2=find(pheta<cut);
        pheta(ibn2) = pheta(ibn2)+2*pi;
        plot(T,pheta)
        xlabel('Time (1 = 4.3 days)')
        ylabel('Thrust Angle (radians)')
        hold on
        
        figure
        J = -.5*jacobi_C(M(:,1:4));
        Tend = T(end);
        plot(T,J)
        title('Jacobi Integral of Motion along Trajectory')
        hold on
        xlabel('T (1 = 4.3 days)')
        ylabel('Jacobi_C')
        L(:,1) = [1.15559740258988      0 0 0]';
        L(:,2) = [0.83702354452395      0 0 0]'; 
        L(:,3) = [-1.00505347015937     0 0 0]';
        L(:,4) = [0.48787143723469   0.86602540378444   0 0]';
        L(:,5) = [0.48787143723469  -0.86602540378444   0 0]';
        JCL=-.5*jacobi_C(L');
        plot([0,Tend],[1,1]*JCL(1),'k') % L_1
        plot([0,Tend],[1,1]*JCL(2),'m') % L_2
        plot([0,Tend],[1,1]*JCL(3),'g') % L_3
        plot([0,Tend],[1,1]*JCL(4),'r') % L_4
        %plot([0,Tend],[1,1]*JCL(5),'m') % L_5
        legend('Trajectory J_C', 'L_1','L_2','L_3','L_4 & L_5')
        figure
        plot3(M(:,1),M(:,2),J)
        title('Trajectory (x(t), y(t), J(X(t))')
        xlabel('x'),ylabel('y'),zlabel('-.5*Jacobi Constant of Motion')


    end

   % hold off
end

global FINPT
FINPT = Mf;
d=d';

end

function S_dot=state_dot_OC(t,S)

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
global TH
Thrust = TH;


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
%     pheta=datan2(-lam4,-lam3)
      l34n = sqrt(lam3*lam3+lam4*lam4);
      cth = -lam3/l34n;
      sth = -lam4/l34n;


% Formulation of State ODEs for x, y xdot, ydot:
      um1or13 = um1/r13;
      uor23 = u/r23;
      S_dot(1) = vx;
      S_dot(2) = vy;
      S_dot(3) = 2*vy + x - um1or13*xpu - uor23*xpum1 + Thrust*cth;
      S_dot(4) =-2*vx + y - y*(um1or13   + uor23)     + Thrust*sth;

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



end
% ___________________________________________________
function [VALUE,ISTERMINAL,DIRECTION] = events(T,X)

global RBOUNDSQ

%VALUE =RBOUNDSQ - X(1:4)'*X(1:4)  ;
VALUE =RBOUNDSQ - X'*X  ;
ISTERMINAL = 1;
DIRECTION = 0;

end
