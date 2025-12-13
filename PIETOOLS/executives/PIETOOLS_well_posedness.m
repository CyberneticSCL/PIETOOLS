%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_wellposedness.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a well-posedness test for a 4-PIE System defined
% by the 2 4-PI operator representation
%   d/dt Top x(t) = Aop x(t)
%
% In particular, we test whether there exists Pop>=eppos*I, Rop>=I s.t.
%   Top'*Pop*Aop + Aop'*Pop*Top <= 2*omega*Top'*Pop*Top     (dissipativity)
%   (Top-Aop)*Rop*(Top-Aop)' >= eppos*I                     (surjectivity)
% Then the pair (Top,Aop) generate a C0-semigroup e^{tA} that satisfies
%   ||e^{tA}|| <= M*e^{omega*t}
% for some M>=1, and where A=Aop*Top^{-1}.
% Note that omega can be positive, and is specified as
%   settings.epneg = -omega;
%
% INPUT: 
% PIE - A pie_struct class object with the above listed PI operators as fields
% settings - An lpisettings() structure with relevant optimization parameters defined
% 
% OUTPUT:
% prog - a solved lpi program structure
% P,R - opvar objects representing operators certifying well-posedness
% omega - if feasible, rate such that ||e^{tA}|| <= M*e^{omega*t}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2025 PIETOOLS Team
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ - 12/11/2025: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prog, P, R, omega] = PIETOOLS_well_posedness(PIE, settings)

% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Pass to the 2D executive if necessary.
if PIE.dim==2
    % Call the 2D version of the executive.
    error("Well-posedness analysis of 2D PDEs is currently not supported.")
    if nargin==1
        [prog, P, R, omega] = PIETOOLS_well_posedness_2D(PIE);
    else
        [prog, P, R, omega] = PIETOOLS_well_posedness_2D(PIE,settings);
    end
    return
end
% Extract PIE operators necessary for the executive.
Top = PIE.T;        Aop = PIE.A;

% get settings information
if nargin<2
    settings_PIETOOLS_heavy;
    settings.sos_opts.simplify = 1; % Use psimplify
    settings.eppos = 1e-2;      % Positivity of Lyapunov Function with respect to real-valued states
    settings.eppos2 = 1*1e-2;   % Positivity of Lyapunov Function with respect to spatially distributed states
    settings.epneg = 0;         % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
end  

dd1 = settings.dd1;
dd12 = settings.dd12;
sos_opts = settings.sos_opts;
options1 = settings.options1;
options12 = settings.options12;
override1 = settings.override1;
eppos = settings.eppos;
omega = -settings.epneg;             % omega is allowed to be positive!
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


fprintf('\n --- Executing Well-Posedness Test --- \n')
% Declare an LPI program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);    % Initialize the program structure


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop and Rop
disp('- Declaring the operator decision variables...');

[prog, P1op] = poslpivar(prog, Top.dim,dd1,options1);
[prog, R1op] = poslpivar(prog, Top.dim,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim,dd12,options12);
    Pop=P1op + P2op;
    [prog, R2op] = poslpivar(prog, Top.dim,dd12,options12);
    Rop=R1op + R2op;
else
    Pop = P1op;
    Rop = R1op;
end

% Enforce strict positivity of the operators
Imat = blkdiag(eppos*eye(Pop.dim(1,:)),eppos2*eye(Pop.dim(2,:)));
Iop1 = mat2opvar(Imat, Pop.dim(:,2), PIE.vars, PIE.dom);
Pop = Pop + Iop1;       % Pop >= eppos*I
Iop2 = mat2opvar(eye(size(Rop)), Rop.dim(:,2), PIE.vars, PIE.dom);
Rop = Rop + Iop2;       % Rop > =I

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the inequality constraints
%
% Dissipativity:     T* P A + A* P T <= 2 omega T* P T
% Surjectvity:       (T-A) R (T*-A*) >= eppos*I

disp('- Constructing the operator inequality constraints...');

Dop1 = Top'*Pop*Aop + Aop'*Pop*Top -2*omega*Top'*Pop*Top;
Dop2 = (Top-Aop)*Rop*(Top-Aop)' - Iop1;
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the operator inequality constraints...');

if sosineq_on
    disp('  - Using lpi_ineq...');
    prog = lpi_ineq(prog,-Dop1,opts);
    prog = lpi_ineq(prog,Dop2,opts);
else
    disp('  - Using an Equality constraint...');
    
    [prog, De1op1] = poslpivar(prog,Dop1.dim,dd2,options2);
    [prog, De1op2] = poslpivar(prog,Dop2.dim,dd2,options2);
    
    if override2~=1
        [prog, De2op1] = poslpivar(prog,Dop1.dim,dd3,options3);
        Deop1 = De1op1+De2op1;
        [prog, De2op2] = poslpivar(prog,Dop2.dim,dd3,options3);
        Deop2 = De1op2+De2op2;
    else
        Deop1 = De1op1;
        Deop2 = De1op2;
    end
    prog = lpi_eq(prog,Dop1+Deop1,'symmetric'); %Dop=-Deop
    prog = lpi_eq(prog,Dop2-Deop2,'symmetric');
end

disp('- Solving the LPI using the specified SDP solver...');
%solving the LPI program
prog = lpisolve(prog,sos_opts); 

% Conclusion:
P = getsol_lpivar(prog,Pop);
R = getsol_lpivar(prog,Rop);
end
