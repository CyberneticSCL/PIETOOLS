%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_stability_dual.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function executes a stability analysis for the dual PIE System defined
% by the 4-PI operator representation
% Top^* \dot x(t)=Aop^* x(t)
%
% NOTE: Stability of the dual system is known to be equivalent to stability of
% the primal system.
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP,SS - 10_01_2020
%  MP - 5_30_2021; changed to new PIE data structure
% SS - 6/1/2021; changed to function, added settings input
% DJ - 06/02/2021; incorporate sosineq_on option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prog, P] = PIETOOLS_stability_dual(PIE, settings)

if PIE.dim==2
    % Call the 2D version of the executive.
    if nargin==1
        [prog, P] = PIETOOLS_stability_dual_2D(PIE);
    else
        [prog, P] = PIETOOLS_stability_dual_2D(PIE,settings);
    end
    return
end

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


% Dumping relevant 4-PI operators to the workspace -MP, 5/2021
Aop=PIE.A;
Top=PIE.T;
fprintf('\n --- Executing Dual Stability Test --- \n')

% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varlist = [Aop.var1; Aop.var2];  % retrieving the names of the independent pvars from Aop (typically s and th)
prog = sosprogram(varlist);      % Initialize the program structure
X=Aop.I;                         % retrieve the domain from Aop
nx1=Aop.dim(1,1);                % retrieve the number of ODE states from Aop
nx2=Aop.dim(2,1);                % retrieve the number of distributed states from Aop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the Lyapunov 
% function candidate
disp('- Parameterizing Positive Lyapunov Operator using specified options...');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the Lyapunov Inequality
%
% i.e. - Assemble the big operator
% Pheq = [ A*P*T'+T*P*A']

disp('- Constructing the Negativity Constraint...');

Dop = [Top*Pop*(Aop')+Aop*Pop*(Top')+epneg*Top*Top']; 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Negativity Constraint...');

% derivative negativity
% constraints
if sosineq_on
    disp('  - Using lpi_ineq...');
    prog = lpi_ineq(prog,-Dop,opts);
else
    disp('  - Using an Equality constraint...');

    [prog, De1op] = poslpivar(prog, [nx1, nx2],X,dd2,options2);
    
    if override2~=1
        [prog, De2op] = poslpivar(prog,[nx1, nx2],X, dd3,options3);
        Deop=De1op+De2op;
    else
        Deop=De1op;
    end
    prog = lpi_eq(prog,Deop+Dop); %Dop=-Deop
end

disp('- Solving the LPI using the specified SDP solver...');
%solving the sos program
prog = sossolve(prog,sos_opts); 
% Conclusion:
P = getsol_lpivar(prog,Pop);
if isfield(settings.sos_opts,'solver')&&strcmp(settings.sos_opts.solver,'sdpt3')
    if exist('prog', 'var')
        if ~(prog.solinfo.info.pinf||prog.solinfo.info.dinf)
            disp('The System of equations was successfully solved.')
        elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf
            disp('The System of equations was not solved.')
        else
            disp('Unable to definitively determine feasibility. Numerical errors dominating or at the limit of stability.')
        end
    end
elseif isfield(settings.sos_opts,'solver')&&strcmp(settings.sos_opts.solver,'sdpnalplus')
    if exist('prog', 'var')
        if ~(prog.solinfo.info.pinf||prog.solinfo.info.dinf)
            disp('The System of equations was successfully solved.')
        elseif ~(prog.solinfo.info.pinf||prog.solinfo.info.dinf) && prog.solinfo.info.numerr
            disp('The System of equations was successfully solved. However, Double-check the precision.')
        elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || prog.solinfo.info.numerr
            disp('The System of equations was not solved.')
        else
            disp('Unable to definitively determine feasibility. Numerical errors dominating or at the limit of stability.')
        end
    end
elseif exist('prog', 'var')
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
    disp('ODE-PDE converted to PIE. No problem solved because executive file was not selected');
end
end
