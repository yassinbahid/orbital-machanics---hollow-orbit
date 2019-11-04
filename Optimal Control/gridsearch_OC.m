function [Xbest, dbest, FinPt,DDS] = gridsearch_OC(maxerr,Ngrid,c1_int,c2_int,k1_int,k2_int)
% Program by Bruce Lundberg, PHD, and Yassin Bahid
% Program gridsearch_5D_adj.m: This program will utilize the shooting
% method defined for our full system in the file shootingF_OC.m. The
% program will attempt to find a value of X for which the error d is below a
% certain user-defined minimum defined in maxerr.

% Ngrid is the number of investigation intervals per parameter.

% ATTENTION ! Current code demands that the user declare PLOTON global and
% enter PLOTON=1, if output plots are desired.
% Also, don't start the time interval [t_begin,t_end] with t_begin=0. (this
% will cause an error in the ode solver)

% Sample input:
% Init_OCP;
% [Xbest, dbest] = gridsearch_OC(1,3,[-.2,.2],[-.2,.2],[-.2,.2],[-.2,.2],[.1,1])



global PLOTON
plotons = PLOTON;
PLOTON = 0;

solutionsX=[];

%Initial Values
c1 = c1_int(1);
c2 = c2_int(1);
k1 = k1_int(1);
k2 = k2_int(1);


% Set Stepsizes

tc1 = (c1_int(2) - c1_int(1));  % Ranges of values. 
tc2 = (c2_int(2) - c2_int(1));
tk1 = (k1_int(2) - k1_int(1));
tk2 = (k2_int(2) - k2_int(1));


if tc1 <= 0
    Ngrid1 = 0;
else
    Ngrid1 = Ngrid;
    tc1 = tc1/Ngrid1;       % dividing up ranges of values into the 
end                         % appropriate number of intervals. 
if tc2 <= 0
    Ngrid2 = 0;
else
    Ngrid2 = Ngrid;
    tc2 = tc2/Ngrid2;
end
if tk1 <= 0
    Ngrid3 = 0;
else
    Ngrid3 = Ngrid;
    tk1 = tk1/Ngrid3;
end
if tk2 <= 0
    Ngrid4 = 0;
else
    Ngrid4 = Ngrid;
    tk2 = tk2/Ngrid4;
end


% Prep plotting DDS vector
NS = (Ngrid1+1)*(Ngrid2+1)*(Ngrid3+1)*(Ngrid4+1);
DDS = zeros([1,NS]);
ls = 0;

dbest = 1.e10;
if isempty(solutionsX)
  
        c1 = c1_int(1);
        for k=1:Ngrid1+1    % if Ngrid1=0, then this loop component will only preform two cycles.
            c2 = c2_int(1);
            for j=1:Ngrid2+1
                k1 = k1_int(1);
                for l=1:Ngrid3+1
                    k2 = k2_int(1);
                    for m=1:Ngrid4+1
                        %                     tf = tf_int(1);
                        %                     for q=1:Ngrid5+1
                        X= [c1,c2,k1,k2]'
                        d = shootingF_OC2(X); dd = norm(d); ls = ls +1;DDS(ls) = dd;
                        %disp([ls,X',dd]);
                        if dd < dbest       % The program will output those values norm error which seem
                            dbest = dd;     % plausible in order to inform the user how close the selected
                            Xbest = X;      % intervals are actually coming to a target solution X.
                            %                             disp('a lower error found:   '),disp(dd)
                            %                             disp('for trial solution :');
                            disp([X',dd]);
                            PLOTON = plotons;
                            if PLOTON ~= 0
                            dd = shootingF_OC2(X);
                            end
                            PLOTON = 0;
                        end
%                         if dd < maxerr      % Here the program will output any solutions X within set tolerance
%                             disp('A solution within specified error tolerances is:');
%                             disp([X',dd]);
%                             %                         end
%                             %                         tf = tf + tt;
%                         end
                        
                        k2 = k2 + tk2;
                    end
                    k1 = k1 + tk1;
                end
                c2 = c2 + tc2;
            end
            c1 = c1 + tc1;
        end
   
end

% % Recompute best Solution with plot
PLOTON = plotons;
dbest = shootingF_OC2(Xbest)'; 
global FINPT
FinPt = FINPT;