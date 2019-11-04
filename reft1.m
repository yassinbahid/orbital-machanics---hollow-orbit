function [rn,vn] = reft1(r,v,mu,t1)
     [a,e,E,i,o,O,nu,tau,A,B] = vec2orbElem(r,v,mu);
     [r1,v1] = ctenewpos(r,v,t1 );
     c = a*(1-e^2);
     f1 = acos((c/norm(r)- 1)/e);
    
%      if e > 1
%         U=@(u) M- e*sinh(u)-u;
%         E = fzero(U, 0);
%         nu = 2*atan(sqrt((e+1)/(e-1))*tanh(E/2));
%      else
%         U=@(u) M-u-e*sin(u);  
%         E = fzero(U, 0);
%         nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
%      end
     r1v = [a*cos(nu)/(1+e*cos(nu));a*sin(nu)/(1+e*cos(nu));0];
     v1v = [-sqrt(mu/a)*sin(nu);-sqrt(mu/a)*(1-cos(nu));0];
     I= [1,0,0; 0, cos(-i),-sin(-i); 0, sin(-i), cos(-i)];
     Om= [cos(-O),-sin(-O),0;sin(-O), cos(-O),0;0,0,1];
     om= [cos(-o),-sin(-o),0;sin(-o), cos(-o),0;0,0,1];
     rn = Om*I*om*r1v;
     vn = Om*I*o*v1v;
     
     