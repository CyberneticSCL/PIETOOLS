%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% executive_PIETOOLS_stability_2D.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a stability analysis for a 2D-PIE System defined
% by the PI operator representation
% Top \dot{u}(t)=Aop u(t)
%
% If any other parts of the PIE are present, these are ignored. Both Top
% and Aop must be properly defined for the script to function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be defined externally:
%
% PIE - PIE data structure. Includes T,A - opvar2ds, typically defined by the conversion script
%
% settings - a matlab structure with following fields are needed, if
% undefined default values are used
% 
% eppos,eppos2,epneg % stricness terms typically defined by the solver script
%
% sos_opts - options for the SOSSOLVER (e.g. sdp solver), typically defined by the solver script
%
% dd1,dd2,dd3,opts,options1,options2,options - accuracy settings, typically defined by the settings script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP, SS, DJ - 28_09_2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prog, Pop, Dop] = executive_PIETOOLS_stability_2D(PIE,settings)

% get settings information
if nargin<2
    settings_PIETOOLS_heavy;
    settings.sos_opts.simplify = 1; % Use psimplify
    settings.eppos = 1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
    settings.eppos2 = 1*1e-6;   % Positivity of Lyapunov Function with respect to spatially distributed states
    settings.epneg = 0*1e-5;    % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
end    
    
dd1 = settings.dd1;
dd12 = settings.dd12;
sos_opts = settings.sos_opts;
options1 = settings.options1;
options12 = settings.options12;
override1 = settings.override1;
eppos = settings.eppos;
epneg = settings.epneg;
eppos2 = settings.eppos2;
ddZ = settings.ddZ;
sosineq_on = settings.sosineq_on;
if sosineq_on
    opts = settings.opts;
else
    override2 = settings.override2;
    options2 = settings.options2;
    options3 = settings.options3;
    dd2 = settings.dd2;
    dd3 = settings.dd3;
end

% Convert degrees to format for poslpivar_2d
dd1 = struct('dx',{dd1});
dd12 = struct('dx',{dd12});
dd2 = struct('dx',{dd2});
dd3 = struct('dx',{dd3});

% Get rid of options which cause errors in poslpivar_2d
options1.exclude = zeros(1,16);  options1.diag = 0;  options1.sep = zeros(1,6);
options12.exclude = zeros(1,16);  options12.diag = 0;  options12.sep = zeros(1,6);
options2.exclude = zeros(1,16);  options2.diag = 0;  options2.sep = zeros(1,6);
options3.exclude = zeros(1,16);  options3.diag = 0;  options3.sep = zeros(1,6);

disp('Executing Primal Stability Test')
% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Aop = PIE.A;    Top = PIE.T;
varlist = [Aop.var1; Aop.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
dom = Aop.I;                       % retrieve the domain from Aop
no1 = Aop.dim(1,1);                % retrieve the number of ODE states from Aop
np2 = Aop.dim(4,1);                % retrieve the number of distributed states from Aop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the Lyapunov 
% function candidate
disp('Parameterizing Positive Lyapunov Operator using specified options');

[prog, P1op] = poslpivar_2d(prog, [no1,0,0,np2],dom,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar_2d(prog, [no1,0,0,np2],dom,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% Enforce strict positivity on the operator
Pop.R00 = Pop.R00 + eppos*eye(no1);
Pop.R22{1,1} = Pop.R22{1,1} + eppos2*eye(np2);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the Lyapunov Inequality
%
% i.e. - Assemble the big operator
% Pheq = [ A'*P*T+T'*P*A]

disp('Constructing the Negativity Constraint');

if epneg==0
    Dop = Top'*Pop*Aop + Aop'*Pop*Top; 
else
    Dop = Top'*Pop*Aop + Aop'*Pop*Top + epneg*(Top'*Top); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('Enforcing the Negativity Constraint');

if sosineq_on
    disp('Using lpi_ineq');
    prog = lpi_ineq_2d(prog,-Dop,opts);
else
    disp('Using an equality constraint');
    
    [prog, De1op] = poslpivar_2d(prog, [no1,0,0,np2],dom,dd2,options2);
    if override2~=1
        [prog, De2op] = poslpivar_2d(prog,[no1,0,0,np2],dom,dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    prog = lpi_eq_2d(prog,Deop+Dop); %Dop=-Deop
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Solve the problem, and extract the solution
%
disp('Solving the LPI using the specified SDP solver');
prog = sossolve(prog,sos_opts); 

% Conclusion:
if exist('prog', 'var')
    if norm(prog.solinfo.info.feasratio-1)<=.3 && ~prog.solinfo.info.numerr
        disp('The System of equations was successfully solved.')
    elseif norm(prog.solinfo.info.feasratio-1)<=.3 && prog.solinfo.info.numerr
        disp('The System of equations was successfully solved. However, Double-check the precision.')
    elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || norm(prog.solinfo.info.feasratio+1)<=.1
        disp('The System of equations was not solved.')
    else
        disp('Unable to definitively determine feasibility. Numerical errors dominating or at the limit of stability.')
    end
else
    disp('System converted to PIE. No problem solved because executive file was not selected');
end
end