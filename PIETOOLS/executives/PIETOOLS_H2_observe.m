%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_H2_observe.m     PIETOOLS 2024
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
% 
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

function [prog, Lop, gam, P, Z] = PIETOOLS_H2_observe(PIE, settings)

if PIE.dim==2
    error('Optimal Estimation of 2D PIEs is currently not supported.')
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

fprintf('\n --- Executing Search for H_2 Optimal Estimator --- \n')

% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varlist = [PIE.A.var1; PIE.A.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
X=PIE.A.I;                         % retrieve the domain from Aop
nx1=PIE.A.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=PIE.A.dim(2,1);                % retrieve the number of distributed states from Aop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The most common usage of this script is to find the minimum hinf gain bound
% In this case, we define the hinf norm variable which needs to be minimized
dpvar gam;
prog = sosdecvar(prog, gam); %this sets gamma as decision var
%prog = sosineq(prog, gam); %this ensures gamma is lower bounded
prog = sossetobj(prog, gam); %this minimizes gamma, comment for feasibility test
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the storage 
% function candidate and the indefinite operator Zop, which is used to
% contruct the estimator gain
disp('- Declaring Positive Storage Operator variable and indefinite Estimator operator variable using specified options...');

[prog, P1op] = poslpivar(prog, PIE.T.dim(:,1),X,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, PIE.T.dim(:,1),X,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% enforce strict positivity on the operator
Pop.P = Pop.P+eppos*eye(nx1);
Pop.R.R0 = Pop.R.R0+eppos*eye(nx2);  

[prog,Zop] = lpivar(prog,[PIE.T.dim(:,1),PIE.C2.dim(:,1)],X,ddZ);
[prog,Wop] = lpivar(prog,[PIE.B1.dim(:,2),PIE.B1.dim(:,2)],X,ddZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% opvar Iw Iz;
% Iw.dim = [PIE.B1.dim(:,2),PIE.B1.dim(:,2)];
% Iz.dim = [PIE.C1.dim(:,1),PIE.C1.dim(:,1)];
% Iw.P = eye(size(Iw.P)); Iz.P = eye(size(Iz.P));
% Iw.R.R0 = eye(size(Iw.R.R0)); Iz.R.R0 = eye(size(Iz.R.R0));


Dop = [(Pop*PIE.A+Zop*PIE.C2)'*(PIE.T)+PIE.T'*(Pop*PIE.A+Zop*PIE.C2)+PIE.C1'*PIE.C1+epneg*PIE.T'*PIE.T]; 

Dop2 = [Wop -(PIE.B1'*Pop+PIE.D21'*Zop'); -(PIE.B1'*Pop+PIE.D21'*Zop')' Pop];

traceVal = trace(Wop.P);
% ensuring scalar inequality gam>trace
prog = sosineq(prog, gam-traceVal);

disp('- Parameterize the derivative inequality...');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Negativity Constraint...');
    disp('  - Using an Equality constraint...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% The old way, where you have to specify everything
    [prog, De1op] = poslpivar(prog, PIE.T.dim(:,1),X,dd2,options2);
    [prog, De3op] = poslpivar(prog, Dop2.dim(:,1),X,dd2,options2);
   [prog, De2op] = poslpivar(prog,PIE.T.dim(:,1),X, dd3,options3);
   [prog, De4op] = poslpivar(prog, Dop2.dim(:,1),X,dd3,options3);
   Deop=De1op+De2op;
   Deop2=De3op+De4op;
    % derivative negativity
    % constraints
    prog = lpi_eq(prog,Deop+Dop,'symmetric'); %Dop=-Deop
    prog = lpi_eq(prog,Deop2-Dop2,'symmetric'); %Dop2=Deop2


%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = sossolve(prog,sos_opts); 

disp('The closed-loop H-2 norm of the given system is upper bounded by:')
if ~isreal(gam)
    disp(sqrt(double(sosgetsol(prog,gam)))); % check the H2 norm, if solved successfully
else 
    disp(gam);
end

gam = sqrt(double(sosgetsol(prog,gam)));

P = getsol_lpivar(prog,Pop);
Z = getsol_lpivar(prog,Zop);

Lop = getObserver(P,Z);
% end

end


