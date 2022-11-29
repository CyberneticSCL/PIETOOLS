%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_Hinf_control.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function executes a synthesis code for H-infty optimal full-state 
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
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding MP,SS - 10_01_2020
%  MP - 5_30_2021; changed to new PIE data structure
% SS - 6/1/2021; changed to function added settings input
% DJ - 06/02/2021; incorporate sosineq_on option, replacd gamma with gam to
%                   avoid conflict with MATLAB gamma function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
%
% PIE - PIE data structure. Include elements T,A,B1,B2,C1,D11,D12 which are 4-PI operators, typically defined by the conversion script
%
% settings - a matlab structure with following fields are needed, if
% undefined default values are used
%
% sos_opts - options for the SOSSOLVER (e.g. sdp solver), typically defined by the solver script
%
% dd1,dd2,dd3,ddZ,opts,options1,options2,options - accuracy settings, typically defined by the settings script

function [prog, Kop, gam, P, Z] = PIETOOLS_Hinf_control(PIE, settings)

if PIE.dim==2
    error('Optimal control of 2D PIEs is currently not supported.')
end

% get settings information
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

fprintf('\n --- Executing Search for H_infty Optimal Controller --- \n')

% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varlist = [PIE.A.var1; PIE.A.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
X=PIE.A.I;                         % retrieve the domain from Aop
nx1=PIE.A.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=PIE.A.dim(2,1);                % retrieve the number of distributed states from Aop
nw=PIE.B1.dim(1,2);                % retrieve the number of real-valued disturbances
nz=PIE.C1.dim(1,1);                % retrieve the number of real-valued regulated outputs
nu=PIE.B2.dim(1,2);                % retrieve the number of real-valued inputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing to see if there are inputs and disturbances at the boundary
if ~(PIE.Tw==0)
    error('Hinf-dual LPI cannot currently be solved for systems with disturbances at the boundary');
end
if ~(PIE.Tu==0)
    error('Hinf-dual LPI cannot currently be solved for systems with inputs at the boundary');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The most common usage of this script is to find the minimum hinf gain bound
% In this case, we define the hinf norm variable which needs to be minimized
dpvar gam;
prog = sosdecvar(prog, gam); %this sets gamma as decision var
prog = sosineq(prog, gam); %this ensures gamma is lower bounded
prog = sossetobj(prog, gam); %this minimizes gamma, comment for feasibility test
%
% Alternatively, the above 3 commands may be commented and a specific gain
% test specified by defining a specific desired value of gamma. This
% results in a feasibility test instead of an optimization problem.
% gamma = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the storage 
% function candidate and the indefinite operator Zop, which is used to
% contruct the estimator gain
disp('- Declaring Positive Storage Operator variable and indefinite Controller operator variable using specified options...');

[prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, [nx1 ,nx2],X,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos2*eye(nx2);  

[prog,Zop] = lpivar(prog,[nu nx1;0 nx2],X,ddZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the KYP matrix
%
% i.e. - Assemble the big operator
% Pheq = [ -gamma*I  D11       (C1*P+D12*Z)*T'
%          D11'         -gamma*I            B1'
%          T*(C1*P+D12*Z)      B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']

Dop = [-gam*eye(nz)   PIE.D11              (PIE.C1*Pop+PIE.D12*Zop)*(PIE.T');
        PIE.D11'                -gam*eye(nw)  PIE.B1';
        PIE.T*(PIE.C1*Pop+PIE.D12*Zop)'       PIE.B1         (PIE.A*Pop+PIE.B2*Zop)*(PIE.T')+PIE.T*(PIE.A*Pop+PIE.B2*Zop)']; 

disp('- Parameterize the derivative inequality...');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Negativity Constraint...');
if sosineq_on
    disp('  - Using lpi_ineq...');
    prog = lpi_ineq(prog,-Dop,opts);
else
    disp('  - Using an Equality constraint...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% The old way, where you have to specify everything
    [prog, De1op] = poslpivar(prog, [nw+nz+nx1, nx2],X,dd2,options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog,[nw+nz+nx1, nx2],X, dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    % derivative negativity
    % constraints
    prog = lpi_eq(prog,Deop+Dop); %Dop=-Deop
end

%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = sossolve(prog,sos_opts); 

disp('The closed-loop H-infty norm of the given system is upper bounded by:')
if ~isreal(gam)
    disp(double(sosgetsol(prog,gam))); % check the Hinf norm, if the solved successfully
else 
    disp(gam);
end
disp('The feedback control gains have been saved in the PI operator K')

gam = double(sosgetsol(prog,gam));

P = getsol_lpivar(prog,Pop);
Z = getsol_lpivar(prog,Zop);
% err = isequal(P.R.R1,P.R.R2);
% if ~all(err(:))
%     Kop = 'Inverse for opvar can only be calculated for opvar with R1=R2';
% else
Kop = getController(P,Z);
% end

end


