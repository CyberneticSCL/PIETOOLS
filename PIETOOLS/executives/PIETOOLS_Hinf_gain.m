%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_Hinf_gain.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes an H-infty gain analysis for a 4-PIE System defined
% by the 5 4-PI operator representation
% Eop \dot x(t)=Aop  x(t) + B1op  w(t)
%          z(t)=C1op x(t) + D11op w(t)
% INPUT: 
% PIE - A pie_struct class object with the above listed PI operators as fields
% settings - An lpisettings() structure with relevant optimization parameters defined
% 
% OUTPUT:
% prog - a solved sosprogram structure from SOSTOOLS
% gam - Hinf norm of the system
% P - Lyapunov function parameter used in Hinf gain LPI
% 
%
% NOTE: At present, there is an implicit assumption that TB1op=0;
%
% If any other parts of the PIE are present, these are ignored. Top, Aop,
% B1op, C1op, and D11op must be properly defined for the script to function.
%

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
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP,SS - 10_01_2020
%  MP - 5_30_2021; changed to new PIE data structure
% SS - 6/1/2021; changed to function, added settings input
% DJ - 06/02/2021; incorporate sosineq_on option, replacd gamma with gam to
%                   avoid conflict with MATLAB gamma function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prog, P, gam] = PIETOOLS_Hinf_gain(PIE, settings)

if PIE.dim==2
    % Call the 2D version of the executive.
    if nargin==1
        [prog, P, gam] = PIETOOLS_Hinf_gain_2D(PIE);
    else
        [prog, P, gam] = PIETOOLS_Hinf_gain_2D(PIE,settings);
    end
    return
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

% Dumping relevant 4-PI operators to the workspace -MP, 5/2021
Aop=PIE.A;
Top=PIE.T;
B1op=PIE.B1;    TB1op = PIE.Tw;
C1op=PIE.C1;
D11op=PIE.D11;

fprintf('\n --- Searching for Hinf gain bound using primal KYP lemma --- \n')
% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varlist = [Aop.var1; Aop.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
X=Aop.I;                         % retrieve the domain from Aop
nx1=Aop.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=Aop.dim(2,1);                % retrieve the number of distributed states from Aop
nw=B1op.dim(1,2);                % retrieve the number of real-valued disturbances
nz=C1op.dim(1,1);                % retrieve the number of real-valued regulated outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The most common usage of this script is to find the minimum hinf gain bound
% In this case, we define the hinf norm variable which needs to be minimized
dpvar gam;
prog = sosdecvar(prog, gam); %this sets gam = gamma as decision var
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
% function candidate
disp('- Declaring Positive Lyapunov Operator variable using specified options...');

[prog, P1op] = poslpivar(prog, PIE.T.dim(:,1),X,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, PIE.T.dim(:,1),X,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the KYP matrix
%
% i.e. - Assemble the big operator
% Pheq = [ -gamma*I  D'       B'*P*T
%          D         -gamma*I C
%          T'*P*B      C'       T'*P*A+A'*P*T]

% adding adjustment for infinite-dimensional I/O
opvar Iw Iz;
Iw.dim = [PIE.B1.dim(:,2),PIE.B1.dim(:,2)];
Iz.dim = [PIE.C1.dim(:,1),PIE.C1.dim(:,1)];
Iw.P = eye(size(Iw.P)); Iz.P = eye(size(Iz.P));
Iw.R.R0 = eye(size(Iw.R.R0)); Iz.R.R0 = eye(size(Iz.R.R0));


disp('- Constructing the Negativity Constraint...');

Dop = [-gam*Iw+TB1op'*(Pop*B1op)+(Pop*B1op)'*TB1op     D11op'             (Pop*B1op)'*Top+TB1op'*(Pop*Aop);
        D11op                                              -gam*Iz        C1op;
        Top'*(Pop*B1op)+(Pop*Aop)'*TB1op                    C1op'              (Pop*Aop)'*Top+Top'*(Pop*Aop)];
    
    
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
    [prog, De1op] = poslpivar(prog, Iw.dim(:,1)+Iz.dim(:,1)+PIE.T.dim(:,1),X,dd2,options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog,Iw.dim(:,1)+Iz.dim(:,1)+PIE.T.dim(:,1),X, dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    prog = lpi_eq(prog,Deop+Dop); %Dop=-Deop
end


%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = sossolve(prog,sos_opts); 

disp('The H-infty norm of the given system is upper bounded by:')
if ~isreal(gam)
    disp(double(sosgetsol(prog,gam))); % check the Hinf norm, if the solved successfully
else 
    disp(gam);
end
P = getsol_lpivar(prog,Pop);
gam = double(sosgetsol(prog,gam));
end