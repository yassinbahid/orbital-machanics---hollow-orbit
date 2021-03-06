% Script for Parameter Continuation Method in Optimal Control
% To use the script certain global parameters must be established
% e.g. for the specifica case here, by running >>Init_OCP_1_0
% this also establishes an initial guess for the adjoints lam_0
% at the initial parameter value of the trace.
% The Newton solver is set to run silently.
% The code can be copied and modified to set up and run other parameter
% traces.


% Set Solver Options

options1 = mmfsolve('default');
options1.MaxIter = 50;
options1.Scale = 'on';
options1.FunTol = 1.e-9;
options1.Jacobian = 'finite';
%options1.Jacobian = 'broyden';
options1.Display = 'off';

% Set up parameter trace value
N = 110;  %number of trace steps
N = 30;  %number of trace steps
p_start = .3;
p_start = .21;
p_end = .1;
% p_end = .2;

% Trace wrt L_4 halo injection point
p_start = pi/1.67;

p_end = p_start + 2*pi;
P=linspace(p_start,p_end,N+1);
%P = [linspace(p_start,.65,100),linspace(.605,.5,100),linspace(.505,p_end,100)];

% Initial Guess at P(1)pased in here or
%Lf =    % Or read from main workspace

Lf_P= zeros(length(Lf),N+1);  % set up storage for solutions


% Confirm solution Lf at Initial Parameter Value P(1)

% Initial parameter value
pk = P(1);

% New initial orbit radius
% X0 = circ_init(pk,pi,1);
% New Final orbit 
[XF,VR VI] = init_halo(4,.05,1,pk);
[Lfn,Fval,ExitFlag,Iters]=mmfsolve('shootingF_OC2',Lf,options1);

% Display Solution
disp([num2str(Iters),' Iterations to Solution at p(',num2str(1),') = ', num2str(pk)])
disp(Lfn')
PLOTON = 2;
d = shootingF_OC2(Lfn);
pause(.1)
PLOTON = 0;

% Check if Solved
if ExitFlag<=1
    Lf_P(:,1) = Lfn;
else
    disp(['Parameter Trace Interation 1 failure '])
    return
end
warning off

% BEGIN MAIN LOOP OVER OC CASES AT VARIOUS PARAMETER VALUES

for k = 2:N+1

    % Next parameter value
    pk = P(k);
    % New initial orbit radius
    %    X0 = circ_init(pk,4*pi/3,-1);
%     X0 = circ_init(pk,pi,1);
    [XF,VR VI] = init_halo(4,.05,1,pk);
    % Produce next guess
    switch k
        case 2

            Lf_g = Lf_P(:,1);

        case {3, 4, 5}

            for ir = 1:length(Lf)
                ci = polyfit(P(1:(k-1)),Lf_P(ir,1:(k-1)),k-2);
                Lf_g(ir) = polyval(ci,pk);
            end

        otherwise

            for ir = 1:length(Lf)
                ci = polyfit(P((k-4):(k-1)),Lf_P(ir,(k-4):(k-1)),3);
                Lf_g(ir) = polyval(ci,pk);
            end
    end

    % Approximate Solution Parameter Value P(k)

    [Lf_P(:,k),Fval,ExitFlag,Iters]=mmfsolve('shootingF_OC2',Lf_g,options1);

    if ExitFlag<=1
        disp([num2str(Iters),' Iterations to Solution at p(',num2str(k),') = ', num2str(pk)])
        disp(Lf_P(:,k)')
        PLOTON = 2;
        d = shootingF_OC2(Lf_P(:,k));
        pause(.1)
        PLOTON = 0;
    else
        disp(['Parameter Trace Iteration ',num2str(k),' FAILURE: Flag =  ',num2str(ExitFlag)])
        disp('  ')
    end

end


figure
plot(P,Lf_P(5,:)), xlabel('orbit_1 radius'),ylabel('optimal t_f (1 = 4.3 days)')

