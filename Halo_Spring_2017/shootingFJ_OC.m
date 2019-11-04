function J = shootingFJ(X)
%Jacobian for shootingF function for use with newton_sys solver

d = shootingF_OC(X);
T = 0;  %Dummy value and variable
ethresh=1.e-3;
THRESH = ethresh*ones(size(X));
VECTORIZED = 0;

global FAC
if isempty(FAC)
    FAC = [  ];
end

[J,FAC] = numjac('shootingFFT_OC',T,X,d,THRESH,FAC,VECTORIZED);

