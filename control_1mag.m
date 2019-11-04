function [Umags,Uc,X] = control_1mag(t)

N=length(t);
Uc= zeros([4,N]);
X= zeros([4,N]);
Umags=zeros([1,N]);

for k = 1:N;
Uc(:,k) = control_1(t(k));
%Umags(k)=dot(Uc(3:4,k), Uc(3:4,k))^.5;
Umags(k)=dot(Uc(:,k), Uc(:,k))^.5;
end