%Two Body Problem
u=6.67*10^-11
%Barycenter fct
function g = G([x1,y1,m1],[x2,y2,m2])
g=[(x1*m1+x2*m2)/(m1+m2),(y1*m1+y2*m2)/(m1+m2)]
end

% Motion relative to the Center of mass
% r: distance between the two objects
function OrbitG=R([r,m1,m2,vx1,vy1,vx2,vy2])
    % We first find the distance between the distance between G and the 2
    % objects
    G = solve ([m1*r1+m2*r2==(m1+m2)*r, m1*r1==m2*r2],[r1,r2]);
    K1=u*m2^2/(m1+m2)^2;
    L=
    plot
        