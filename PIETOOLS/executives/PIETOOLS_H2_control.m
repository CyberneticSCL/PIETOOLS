%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_H2_control.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function executes a synthesis code for H-2 optimal full-state 
% feedback controller design (w/o control at the boundary) for a 4-PIE 
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
% K - controller gains that stabilize the system has Hinf performance 
% gam - Hinf norm for the obtained controller
% Wo - Lyapunov function parameter that proves stability
% Z - Controller variable used to linearize the Bilinearity in the H2 LPI
% 
% NOTE: The resulting controller has the form
% u(t) = Zop*(Wop)^{-1} x(t)
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
% DB, 01/07/2025: Correct the LPIs and add error for disturbances at the
% boundary
function [prog, Kop, gam, P, Z, W] = PIETOOLS_H2_control(PIE, settings)
%% Extract PIE operators necessary for the executive.
Top = PIE.T;         B1op = PIE.B1;
Aop = PIE.A;        B2op = PIE.B2;
C1op = PIE.C1;      D12op = PIE.D12;
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
%% Inputs at the boundaries 
if ~(PIE.Tu==0)
    error('H2-dual LPI cannot currently be solved for systems with inputs at the boundary');
end
%% Feedthrought terms
if ~(PIE.D11==0)
    error('H2 LPIs requires no feedthrought term');
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
fprintf('\n --- Executing Search for H_2 Optimal Controller --- \n')
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);      % Initialize the program structure
nx1 = PIE.A.dim(1,1);                % retrieve the number of ODE states from Aop
nx2 = PIE.A.dim(2,1);                % retrieve the number of distributed states from Aop
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The most common usage of this script is to find the minimum hinf gain bound
% In this case, we define the hinf norm variable which needs to be minimized
dpvar gam;
prog = lpidecvar(prog, gam); %this sets gamma as decision var
% Alternatively, the above 3 commands may be commented and a specific gain
% test specified by defining a specific desired value of gamma. This
% results in a feasibility test instead of an optimization problem.
% gamma = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: declare the posopvar variables, Wop and Xop, which defines the storage 
% function candidate and an auxiliare operator used to get rid off the nonlinearity on the LPIs.
% Moreover, the indefinite operators Zop is used to contruct the controller gain
disp('- Declaring Positive Storage Operator variables and indefinite Controller operator variable using specified options...');
[prog, P1op] = poslpivar(prog, [nx1;nx2], dd1, options1);
if override1~=1
    [prog, P2op] = poslpivar(prog, [nx1;nx2], dd12, options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end
%% enforce strict positivity on the Lyapunov operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  
%% declare indefinite auxiliary operator Z
dimZ=[B2op.dim(:,2),Top.dim(:,1)];
[prog,Zop] = lpivar(prog, dimZ, ddZ);
%% declare positive auxiliary operator W
dimW=[C1op.dim(:,1),C1op.dim(:,1)];
[prog,Wop] = lpivar(prog, dimW, ddZ);
%% STEP 2: Using the controlability gramian
disp('- Constructing the Inequality Constraints...');
Iw = mat2opvar(eye(size(Bwop,2)), Bwop.dim(:,2), PIE.vars, PIE.dom);

Dneg=[-gam*Iw          B1op'
            B1op     Top*Pop*Aop'+Aop*Pop*Top'+B2op*Zop*Top'+Top*Zop'*B2op'];
Dp12=C1op*Pop+D12op*Zop;
Dpos=[Wop Dp12
            Dp12' Pop];
traceVal = trace(Wop.P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3: Impose Inequality Constraint
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
%%  ensuring scalar inequality gam>trace
prog = lpi_ineq(prog, gam-traceVal);
prog = lpi_ineq(prog, gam); %this ensures gamma is lower bounded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Objective function
% min gam
% s.j. LPIs
prog = lpisetobj(prog, gam); %this minimizes gamma, comment for feasibility test
%% STEP 4: solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = lpisolve(prog,sos_opts); 
%% Evaluating decision variables
P = lpigetsol(prog,Pop);
W = lpigetsol(prog,Wop);
Z =  lpigetsol(prog,Zop);
gam = double(lpigetsol(prog, gam));
%% Display closed-loop H2 norm
disp('The closed-loop H-2 norm of the given system is upper bounded by:')
disp(gam);
%% Controller Reconstruction
Kop = getController(P,Z);
end


