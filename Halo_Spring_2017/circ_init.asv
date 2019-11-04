function X0 = circ_init(r,theta,dir, body)
% Computes initial conditions for a two body circular orbit
% with radius r DU (DU = dimensionless distance units with 1 = primary to
%                   secondary distance)
% and position theta from positive x axis. (clockwise orbit if dir < 0)
% for about the primary (or secondary if body variable is present)


% Reduced Mass Parameter mu = m1/(m1+m2)
global U

% Which central body?
if nargin > 3
    u =  U;  %  secondary mass
    center = 1-U;
else
    u = 1 - U; % primary mass
    center = -U;
end

% Which central body?
if nargin <3 
    dir =  1;  %  ccl revolution
end

omega = 1;  % angular speed of rotating coordinate system
v = sign(dir)*sqrt(u/r) - omega*r; % velocity from vis-viva equation adjusted for rot. coord. syst.

cth = cos(theta);
sth = sin(theta);

X0 = [r*cth + center; r*sth; -v*sth;  v*cth];

% format long g
% format compact
% x0 = X0(1)
% y0 = X0(2)
% vx0 = X0(3)
% vy0 = X0(4)
% format short g


