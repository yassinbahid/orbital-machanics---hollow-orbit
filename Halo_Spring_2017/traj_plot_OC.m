% OC Solution Traj and Steering Angle Plot Script
% Plot Thrust Arc
PLOTON = 1;
d=shootingF_OC2(Lf);

% Plot Coast After Thrust Arc

global FINPT
tfc = 30;
Pf = Traj_simulation(tfc,FINPT(1:4));

% Plot Coast Without Thrust Arc
Pf = Traj_simulation(tfc,X0);

% Replot Thrust Arc and Steering Angle
subplot(2,1,2)
hold off
PLOTON = 2;
d=shootingF_OC2(Lf);