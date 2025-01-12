%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_H2_estimator.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function executes a synthesis code for H-2 optimal observer design (w/o control at the boundary) for a 4-PIE 
% System defined by the 7 4-PI operator representation
% Top \dot x(t)=Aop  x(t) +  B1op w(t) +  B2op u(t)
%          z(t)=C1op x(t) + D11op w(t) + D12op u(t)
%
% INPUT: 
% PIE - A pie_struct class object with the above listed PI operators as fields
% settings - An lpisettings() structure with relevant optimization parameters defined
% 
% OUTPUT:
% prog - a solved sosprogram structure from SOSTOOLS
% L - controller gains that stabilize the system has Hinf performance 
% gam - Hinf norm for the obtained controller
% P - Lyapunov function parameter that proves stability
% Z - Controller variable used to linearize the Bilinearity in the Hinf LPI
% 
% NOTE: The resulting controller has the form
% u(t) = Zop*(Pop)^{-1} x(t)
%
% NOTE: At present, there is an implicit assumption that TB2op=0;
%
% If any other parts of the PIE are present, these are ignored. Top, Aop,
% B1op, B2op, C1op, D11op and D12op must be properly defined for the script to function.
%
% NOTE: Optimizing the Hinf gain bound for the controller synthesis problem
% sometimes results in large errors in sedumi. It is unknown
% how to predict when this will occur. Future releases will address this
% problem. At present, if this occurs, we recommend using bisection on a
% fixed Hinf gain bound or searching for a stabilizing controller.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DB, 2024: Initial coding;
% DJ, 11/30/2024: Update to use LPI programming structure;
% DJ, 01/06/2025: Rename to "estimator";
% DB, 01/07/2025: Correct the LPIs and add error for disturbances at the
% boundary
 
function [prog, Lop, gam, P, Z, W] = PIETOOLS_H2_estimator(PIE, settings)
%% Extract PIE operators necessary for the executive.
Top = PIE.T;         B1op = PIE.B1;
Aop = PIE.A;        C2op = PIE.C2;
C1op = PIE.C1;      D21op = PIE.D21;
Twop=PIE.Tw;
%% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
%%  2D case
if PIE.dim==2
    error('Optimal Estimation of 2D PIEs is currently not supported.')
end
%% Disturbances at the boundaries
if ~(Twop==0)
    error('H2 estimator LPIs cannot currently be solved for systems with disturbances at the boundary');
end
%% Feedthrough terms
if ~(PIE.D11==0)
    error('H2 LPIs requires no feedthrough term');
end
%% get settings information
if nargin<2
    settings_PIETOOLS_light;
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
%% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n --- Executing Search for H_2 Optimal Estimator --- \n')
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);      % Initialize the program structure
nx1=PIE.A.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=PIE.A.dim(2,1);                % retrieve the number of distributed states from Aop
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The most common usage of this script is to find the minimum hinf gain bound
% In this case, we define the hinf norm variable which needs to be minimized
dpvar gam;
prog = lpidecvar(prog, gam); %this sets gamma as decision var
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: declare the posopvar variable, Pop, which defines the storage 
% function candidate and the indefinite operator Zop, which is used to
% contruct the estimator gain
disp('- Declaring Positive Storage Operator variable and indefinite Estimator operator variable using specified options...');
[prog, P1op] = poslpivar(prog, PIE.T.dim(:,1), dd1, options1);
if override1~=1
    [prog, P2op] = poslpivar(prog, PIE.T.dim(:,1), dd12, options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end
%% enforce strict positivity on the Lyapunov operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  
%% declare indefinite auxiliary operator Z
dimZ=[Top.dim(:,1),C2op.dim(:,1)];
[prog,Zop] = lpivar(prog,dimZ, ddZ);
%% declare positive auxiliary operator W
dimW=[B1op.dim(:,2),B1op.dim(:,2)];
[prog,Wop] = lpivar(prog,dimW, ddZ);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Using the observability gramian
disp('- Constructing the Inequality Constraints...');
Iz = mat2opvar(eye(size(C1op,1)), C1op.dim(:,1), PIE.vars, PIE.dom);

Dneg=[-gam*Iz          C1op
            C1op' Top'*Pop*Aop+Aop'*Pop*Top+Top'*Zop*C2op+C2op'*Zop'*Top];
D12=B1op'*Pop+D21op'*Zop';
Dpos=[Wop -D12
            -D12' Pop];
traceVal = trace(Wop.P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Inequalities Constraints...');
if sosineq_on
    disp('  - Using lpi_ineq...');
    prog = lpi_ineq(prog,-Dneg,opts);
    prog = lpi_ineq(prog,Dpos,opts);
else
    disp('  - Using an Equality constraint...');
    [prog, De1op] = poslpivar(prog, Dneg.dim, dd2, options2);
    [prog, De3op] = poslpivar(prog, Dpos.dim, dd2, options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog, Dneg.dim, dd3, options3);
        Deop=De1op+De2op;
         [prog, De4op] = poslpivar(prog, Dpos.dim, dd3, options3);
        Deopp=De3op+De4op;      
    else
        Deop=De1op;
        Deopp=De3op;
    end
    prog = lpi_eq(prog,Deop+Dneg,'symmetric'); %Dneg=-Deop
    prog = lpi_eq(prog,Deopp-Dpos,'symmetric'); %Dpos=Deopp
end
%% ensuring scalar inequality gam>trace
prog = lpi_ineq(prog, gam-traceVal);
prog = lpi_ineq(prog, gam); %this ensures gamma is lower bounded
%% Objective function
% min gam
% sj LPIs
prog = lpisetobj(prog, gam); %this minimizes gamma, comment for feasibility test
%% STEP 4: solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = lpisolve(prog,sos_opts); 
%% Evaluating decision variables
P = lpigetsol(prog,Pop);
W = lpigetsol(prog,Wop);
Z = lpigetsol(prog,Zop);
gam = double(lpigetsol(prog,gam));
%% Display closed-loop H2 norm
disp('The closed-loop H-2 norm of the given system is upper bounded by:')
disp(gam);
%% Observer Reconstruction
Lop = getObserver(P,Z);
end


