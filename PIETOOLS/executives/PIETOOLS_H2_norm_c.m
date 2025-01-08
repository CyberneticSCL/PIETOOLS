function [prog,W,gam,R,Q]= PIETOOLS_H2_norm_c(PIE, settings,options)
% This function solves the minimization problem to obtain the H2 norm of a linear distributed parameter
% system using the controlability gramian approach and PIE framework. For
% the feasibility test an additional input with a assigned value to the
% norm is required.
% inputs: (mandatory)
%   PIE : PIE structure of the corresponding system
%   settings : options related to PIETOOLS
% outputs: 
%   h2 = computed H2 norm for SDP problem or optional input value for feasibility tests;
%   W = controlability gramian, a positive PI operator;
%   prog = sum of squares program information;
%   R = positive PI operator parameterizing the non-coercive Lyapunov 
%       functional V(v)=<v,R*v>;
%   Q = indefinite PI operator such that R=PIE.T*Q=Q'*PIE.T';
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, D. Jagt at djagt@asu.edu, or D.
% Braghini d166353@dac.unicamp.br.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% DJ, 10/19/2024: Update to use new LPI programming structure;
% DJ, 01/06/2025: Update to non-coercive version;
% DJ, 01/07/2025: Correction to output gam, should not take square root;


% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Pass to the 2D executive if necessary.
if PIE.dim==2
    if nargin==1
        [prog, W, gam] = PIETOOLS_H2_norm_2D_c(PIE);
    elseif nargin==2
        [prog, W, gam] = PIETOOLS_H2_norm_2D_c(PIE, settings);
    else
        [prog, W, gam] = PIETOOLS_H2_norm_2D_c(PIE, settings, options);
    end
    return
end
% Extract PIE operators necessary for the executive.
Top = PIE.T;        Twop = PIE.Tw;
Aop = PIE.A;        B1op = PIE.B1;
Czop = PIE.C1;

% Make sure thera are no disturbances at the boundary.
if ~(Twop==0)
    error('H2 norm LPI cannot currently be solved for systems with disturbances at the boundary');
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

fprintf('\n --- Searching for H2 norm bound using the controlability gramian --- \n')
% Declare an SOS program and initialize domain and opvar spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);      % Initialize the program structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if the user wants to calculate the norm, define objective function
% otherwise, the norm is input to the function.
if nargin<=2 || ~isfield(options,'h2')
    dpvar gam;
    prog = lpidecvar(prog, gam); % set gam = gamma as decision variable
    prog = lpi_ineq(prog, gam);  % enforce gamma>=0
    prog = lpisetobj(prog, gam); % set gamma as objective function to minimize
else
    gam = options.h2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Rop, which defines the storage 
% function candidate V(v)=<v,Rop*v>=<Top'*v,Wop*Top'*v>
disp('- Declaring Gramian using specified options...');

[prog, R1op] = poslpivar(prog, Top.dim, dd1, options1);

if override1~=1
    [prog, P2op] = poslpivar(prog, Top.dim, dd12, options12);
    Rop=R1op+P2op;
else
    Rop=R1op;
end

% Also declare an indefinite operator Qop=Wop*Top' so that Rop = Top*Qop.   % DJ, 01/06/2025
Qdeg = get_lpivar_degs(Rop,Top);
[prog, Qop] = lpivar(prog,Top.dim,Qdeg);
prog = lpi_eq(prog, Top*Qop-Rop);

% Finally, declare a positive operator (matrix) Wm representing the
% controllability Gramian
[prog, Wm] = poslpivar(prog,Czop.dim(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Using the controlability gramian

disp('- Constructing the Negativity Constraint...');
Iw = mat2opvar(eye(size(B1op,2)), B1op.dim(:,2), PIE.vars, PIE.dom);

Dneg = [-gam*Iw   B1op'
        B1op      Qop'*Aop'+Aop*Qop];
Dpos = [Wm           Czop*Qop
        Qop'*Czop'   Rop];
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
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

% ensuring scalar inequality gam>trace
traceVal = trace(Wm.P);
prog = lpi_ineq(prog, gam-traceVal);

%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = lpisolve(prog,sos_opts); 

W = lpigetsol(prog,Wm);
R = lpigetsol(prog,Rop);
Q = lpigetsol(prog,Qop);

if nargin<=2 || ~isfield(options,'h2')
        gam = double(lpigetsol(prog,gam));                                  % DJ, 01/07/2025
        disp('The H2 norm of the given system is upper bounded by:')
         disp(gam);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%