%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_Hinf_gain.m     PIETOOLS 2024
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
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP,SS - 10_01_2020
% MP - 05/30/2021: changed to new PIE data structure;
% SS - 06/01/2021: changed to function, added settings input;
% DJ - 06/02/2021: incorporate sosineq_on option, replacd gamma with gam to
%                   avoid conflict with MATLAB gamma function;
% DJ - 10/19/2024: Update to use new LPI programming structure;
% DJ - 01/06/2024: Update to non-coercive version;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prog, R, gam] = PIETOOLS_Hinf_gain(PIE, settings)

% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Pass to the 2D executive if necessary.
if PIE.dim==2
    % Call the 2D version of the executive.
    if nargin==1
        [prog, R, gam] = PIETOOLS_Hinf_gain_2D(PIE);
    else
        [prog, R, gam] = PIETOOLS_Hinf_gain_2D(PIE,settings);
    end
    return
end
% Extract PIE operators necessary for the executive.
Top = PIE.T;    Twop = PIE.Tw;
Aop = PIE.A;    Bwop = PIE.Bw;    
Czop = PIE.Cz;  Dzwop = PIE.Dzw;

% Don't support boundary disturbances
if ~(Twop==0)
        fprintf('\n --- Non-coercive LPIs cannot currently be solved for systems with disturbances at the boundary. Calling the coercive version.---\n');
     if nargin==1
         [prog, R, gam] = PIETOOLS_Hinf_gain_coercive(PIE);
    else
         [prog, R, gam] = PIETOOLS_Hinf_gain_coercive(PIE,settings);
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
eppos2 = settings.eppos2;
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



fprintf('\n --- Searching for Hinf gain bound using primal KYP lemma --- \n')
% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);      % Initialize the program structure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The most common usage of this script is to find the minimum hinf gain bound
% In this case, we define the hinf norm variable which needs to be minimized
dpvar gam;
prog = lpidecvar(prog, gam); % set gam = gamma as decision variable
prog = lpi_ineq(prog, gam);  % enforce gamma>=0
prog = lpisetobj(prog, gam); % set gamma as objective function to minimize
%
% Alternatively, the above 3 commands may be commented and a specific gain
% test specified by defining a specific desired value of gamma. This
% results in a feasibility test instead of an optimization problem.
% gamma = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Rop, which defines the storage 
% function candidate V(v)=<v,Rop*v>=<Top*v,Pop*Top*v>
disp('- Declaring Positive Lyapunov Operator variable using specified options...');

[prog, R1op] = poslpivar(prog, Top.dim,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim,dd12,options12);
    Rop=R1op+P2op;
else
    Rop=R1op;
end

% Also declare an indefinite operator Qop=Pop*Top so that Rop = Top'*Qop.   % DJ, 01/06/2025
Qdeg = get_lpivar_degs(Rop,Top);
[prog, Qop] = lpivar(prog,Top.dim,Qdeg);
prog = lpi_eq(prog, Top'*Qop-Rop);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the KYP matrix
%
% i.e. - Assemble the big operator
% Pheq = [ -gamma*I  D'       B'*P*T
%          D         -gamma*I C
%          T'*P*B      C'       T'*P*A+A'*P*T]

disp('- Constructing the Negativity Constraint...');

Iw = mat2opvar(eye(size(Bwop,2)), Bwop.dim(:,2), PIE.vars, PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)), Czop.dim(:,1), PIE.vars, PIE.dom);

Dop = [-gam*Iw,       Dzwop',     Bwop'*Qop;
        Dzwop,        -gam*Iz,    Czop;
        Qop'*Bwop,    Czop'       Aop'*Qop+Qop'*Aop];
    
    
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
    [prog, De1op] = poslpivar(prog, Dop.dim, dd2, options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog, Dop.dim, dd3, options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    prog = lpi_eq(prog,Deop+Dop,'symmetric'); %Dop=-Deop
end


%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = lpisolve(prog,sos_opts); 

disp('The H-infty norm of the given system is upper bounded by:')
if ~isreal(gam)
    disp(double(lpigetsol(prog,gam))); % check the Hinf norm, if solved successfully
else 
    disp(gam);
end
R = lpigetsol(prog,Rop);
gam = double(lpigetsol(prog,gam));
end