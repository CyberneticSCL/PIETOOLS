%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_Hinf_estimator.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script executes an H-infty gain analysis for a 4-PIE System defined
% by the 7 4-PI operator representation
% Top \dot x(t)=Aop  x(t) + B1op  w(t)
%          z(t)=C1op x(t) + D11op w(t)
%          y(t)=C2op x(t) + D21op w(t)
%
% INPUT: 
% PIE - A pie_struct class object with the above listed PI operators as fields
% settings - An lpisettings() structure with relevant optimization parameters defined
% 
% OUTPUT:
% prog - a solved sosprogram structure from SOSTOOLS
% L - observer gains that stabilize the system has Hinf performance 
% gam - Hinf norm for the obtained observer
% P - Lyapunov function parameter that proves stability of the observer
%     error
% Z - Observer variable used to linearize the Bilinearity in the Hinf LPI
% 
% NOTE: The resulting estimator has the form
% Top \dot \hat x(t)=Aop  \hat x(t) + (Pop)^{-1}Zop*(C2op \hat x(t)-y(t))
%
% NOTE: At present, there is an implicit assumption that TB1op=0;
%
% If any other parts of the PIE are present, these are ignored. Top, Aop,
% B1op, C1op, C2op, D11op, and D21op must be properly defined for the script to function.
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MP,SS - 10_01_2020
% MP - 05/30/2021: changed to new PIE data structure;
% SS - 06/01/2021: changed to function, added settings input;
% DJ - 06/02/2021: incorporate sosineq_on option, replacd gamma with gam to
%                   avoid conflict with MATLAB gamma function;
% DJ - 08/01/2024: add call to 2D version;
% DJ - 10/19/2024: Update to use new LPI programming structure;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prog, Lop, gam, P, Z] = PIETOOLS_Hinf_estimator(PIE, settings,gain)

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
        [prog, Lop, gam, P, Z] = PIETOOLS_Hinf_estimator_2D(PIE);
    elseif nargin==2
        [prog, Lop, gam, P, Z] = PIETOOLS_Hinf_estimator_2D(PIE,settings);
    elseif nargin==3
        [prog, Lop, gam, P, Z] = PIETOOLS_Hinf_estimator_2D(PIE, settings, gain);
    end
    return
end
% Extract PIE operators necessary for the executive.
Top = PIE.T;        Twop = PIE.Tw;
Aop = PIE.A;        Bwop = PIE.B1;    
Czop = PIE.C1;      Dzwop = PIE.D11;
Cyop = PIE.C2;      Dywop = PIE.D21;

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

fprintf('\n --- Executing Search for H_infty Optimal Estimator --- \n')



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
% STEP 1: declare the posopvar variable, Pop, which defines the storage 
% function candidate and the indefinite operator Zop, which is used to
% contruct the estimator gain
disp('- Declaring Positive Storage Operator variable and indefinite Observer operator variable using specified options...');

[prog, P1op] = poslpivar(prog, Top.dim, dd1, options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim, dd12, options12);
    Pop=P1op+P2op;
else
    Pop=P1op;
end

% Enforce strict positivity of the operator
Imat = blkdiag(eppos*eye(Pop.dim(1,:)),eppos2*eye(Pop.dim(2,:)));
Pop = Pop + mat2opvar(Imat, Pop.dim(:,2), PIE.vars, PIE.dom);

[prog,Zop] = lpivar(prog,Cyop.dim(:,[2,1]),ddZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Define the KYP matrix
%
% i.e. - Assemble the big operator
% Pheq = [ -gamma*I  -D11'       -(P*B1+Z*D21)'*T
%          -D11         -gamma*I            C1
%          -T'*(P*B1+Z*D21)      C1'       T'*(P*A+Z*C2)+(P*A+Z*C2)'*T]

disp('- Constructing the Negativity Constraint...');

Iw = mat2opvar(eye(size(Bwop,2)), Bwop.dim(:,2), PIE.vars, PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)), Czop.dim(:,1), PIE.vars, PIE.dom);

Dop = [-gam*Iw+Twop'*(Pop*Bwop+Zop*Dywop)+(Pop*Bwop+Zop*Dywop)'*Twop,   -Dzwop',   -(Pop*Bwop+Zop*Dywop)'*Top-Twop'*(Pop*Aop+Zop*Cyop);
       -Dzwop,                                                          -gam*Iz,   Czop;
       -Top'*(Pop*Bwop+Zop*Dywop)-(Pop*Aop+Zop*Cyop)'*Twop,             Czop',     (Pop*Aop+Zop*Cyop)'*Top+Top'*(Pop*Aop+Zop*Cyop)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Negativity Constraint...');
if sosineq_on
    disp('  - Using lpi_ineq');
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

disp('The H-infty gain from disturbance to error in estimated state is upper bounded by:')
if ~isreal(gam)
    disp(double(lpigetsol(prog,gam))); % check the Hinf norm, if the solved successfully
else 
    disp(gam);
end

gam = double(lpigetsol(prog,gam));

P = lpigetsol(prog,Pop);
Z = lpigetsol(prog,Zop);
Lop = getObserver(P,Z);

end