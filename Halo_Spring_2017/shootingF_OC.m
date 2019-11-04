function d = shootingF_OC(X)

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
options1 = odeset('RelTol',TOLR,'AbsTol',TOLA);
[T,M]=ode45('state_dot_OC',tf,M0,options1);
Mf = M(end,:)';

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
    legend('X_0','X(t_f) Attained', 'X_f Target')
    xlabel('x-axis'), ylabel('y-axis'), title('xy-trajectory plot')
    
%     if x_tf-X0(1) > 1  &   y_tf-X0(2) > 1
%         text(x_tf-.1*x_tf, y_tf-.1*y_tf,['(',num2str(x_tf),' ',num2str(y_tf),') ',num2str(d(3)), '  ', num2str(d(4)),'  ' num2str(d(5))])
%         axis equal
%     else 
%         avg_x = abs(mean( x_tf-X0(1)));
%         avg_y = abs(mean( y_tf-X0(2)));
%         text(avg_x,avg_y,['(',num2str(x_tf),' ',num2str(y_tf),') ',num2str(d(3)), '  ', num2str(d(4)),'  ' num2str(d(5))])
%     end
    
    subplot(2,1,2)
    pheta=atan2(-M(:,8),-M(:,7)); % Thrust angle.
    plot(T,pheta)
    xlabel('Time (1 = 4.3 days)')
    ylabel('Thrust Angle (radians)')
    
    hold off 
end

global FINPT
FINPT = Mf;
d=d';