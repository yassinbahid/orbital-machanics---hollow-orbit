function x = newton_sysdalldamp(F,JF,xx0,tol,maxit,damp,newtons)
% F is the system given by an (n x 1) matrix;
% JF is the Jacobian of F given by an (n x n) matrix;
% xx0 is the (n x 1) initial vector, tol ia a tolerance;
% max_it is the maximum number of iterations.
% damp, if present, scales the step vector
% newtons, if 0, uses the J(xx0) for all iterations
% Also computes and prints the final Jacobian
% and its eigenvalues

x0=xx0;
if nargin <6
    damp = 1;
end
if nargin <7
    newtons = 1;
end
global FAC
FAC = [  ];
x = zeros(length(xx0),maxit+1);
x(:,1)=xx0;

fprintf('%i',0),disp([x0'])
iter=1;
while(iter<=maxit)
    % compute current F(x) and the 'Merit Function' dn = ||F(x)||^2
   yn= feval(F,x0);
   dn = yn'*yn;
   if iter == 1
       errf0 = sqrt(dn);      % Record initial distance of F(x) from 0.
   end
   %fprintf('%i',iter);,disp([yn', sqrt(dn)]);%Displays Components of F(x_n)
   if iter == 1 | newtons
       Jn = feval(JF,x0);
   end
   
   % Newton Step
   y = -Jn\yn;
   xn=x0+y;
  
   if damp<1 
   % Test Newton Step
     ynn= feval(F,xn);
     dnn = yn'*yn;
     dampn = 1.;
   end
   while damp<1 & dnn>=dn & dampn>5.e-5
       
       dampn = .2*dampn;
       xn = x0+dampn*y;
       ynn= feval(F,xn);
       dnn = ynn'*ynn;
       
   end
   
   % Compute new distance of F(x) from 0
   if damp < 1
     errf = sqrt(dnn);
   else
     errf = sqrt(dn);
   end
       
 
   x(:,iter+1)=xn;
   err= max(abs(xn-x0));
%    fprintf('%i',iter);,disp([xn', sqrt(dn)]);
   fprintf('%3i   %4g  %6g ',iter,xn', sqrt(dn));
   fprintf('\n');
   if err<=tol && errf<=tol*max(errf0,.001)
      fprintf('  Newton''s method converges after %3.0f iterations to \n',iter);
      %    x
	  x = x(:,1:(iter+1));
      J = feval(JF,xn);
%       fprintf('\n The Jacobian at x: \n ')
%       disp(J)
%       fprintf('\n Eigenvalues of the Jacobian: ')
%       fprintf(' %g   ',eig(J)')
%       fprintf('\n \n ')
      return;
   else
      x0=xn;
   end
   iter=iter+1;
end
disp('Newton''s method does not converge')
x(:,iter+1)=xn;
