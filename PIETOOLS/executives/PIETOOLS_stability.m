%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_stability.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes a stability analysis for a 4-PIE System defined
% by the 2 4-PI operator representation
% Top \dot x(t)=Aop x(t)
%
% If any other parts of the PIE are present, these are ignored. Both Top
% and Aop must be properly defined for the script to function.
%
% INPUT: 
% PIE - A pie_struct class object with the above listed PI operators as fields
% settings - An lpisettings() structure with relevant optimization parameters defined
% 
% OUTPUT:
% prog - a solved sosprogram structure from SOSTOOLS
% P - Lyapunov function parameter that proves stability

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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP,SS - 10_01_2020
% MP - 05/30/2021: changed to new PIE data structure;
% SS - 06/01/2021: changed to function, added settings input;
% DJ - 06/02/2021: incorporate sosineq_on option;
% DJ - 10/19/2024: Update to use new LPI programming structure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prog, P] = PIETOOLS_stability(PIE, settings)

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
        [prog, P] = PIETOOLS_stability_2D(PIE);
    else
        [prog, P] = PIETOOLS_stability_2D(PIE,settings);
    end
    return
end
% Extract PIE operators necessary for the executive.
Top = PIE.T;        Aop = PIE.A;

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


fprintf('\n --- Executing Primal Stability Test --- \n')
% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);    % Initialize the program structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the Lyapunov 
% function candidate
disp('- Parameterizing Positive Lyapunov Operator using specified options...');

[prog, P1op] = poslpivar(prog, Top.dim,dd1,options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim,dd12,options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% enforce strict positivity of the operator
Imat = blkdiag(eppos*eye(Pop.dim(1,:)),eppos2*eye(Pop.dim(2,:)));
Pop = Pop + mat2opvar(Imat, Pop.dim(:,2), PIE.vars, PIE.dom);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the Lyapunov Inequality
%
% i.e. - Assemble the big operator
% Pheq = [ A'*P*T+T'*P*A]

disp('- Constructing the Negativity Constraint...');

Dop = Top'*Pop*Aop + Aop'*Pop*Top +epneg*Top'*Pop*Top; 
    


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
    
    [prog, De1op] = poslpivar(prog,Dop.dim,dd2,options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog,Dop.dim,dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    prog = lpi_eq(prog,Dop+Deop,'symmetric'); %Dop=-Deop
end

disp('- Solving the LPI using the specified SDP solver...');
%solving the sos program
prog = lpisolve(prog,sos_opts); 

% Conclusion:
P = getsol_lpivar(prog,Pop);
end
function prog = remove_dup(prog)
for i=1:length(prog.expr.At)
    At = prog.expr.At{i};
    b = prog.expr.b{i};
    C = [b';At];
    C = unique(C','rows','stable');
    C = C';
    prog.expr.b{i} = C(1,:)';
    prog.expr.At{i} = C(2:end,:);
    prog.expr.Z{i} = unique(prog.expr.Z{i},'rows');
end
end
