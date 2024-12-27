function [prog,W, gam, R,Q]= PIETOOLS_H2_norm_o(PIE, settings,options)
% This function solves a minimization problem to obtain the H2 norm of a linear distributed parameter
% system using the observability gramian approach and PIEs framework. For
% the feasibility test an additional input with a assigned value to the
% norm is required.
% inputs: (mandatory)
%   PIE : PIE structure of the corresponding system
%   settings : options related to PIETOOLS
% outputs:
%   gam = computed H2 norm for SDP problem or optional input value for feasibility tests;
%   Wo= Observability gramian, a positive PI operator;
%  prog= sum of squares program information;
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, D. Jagt at djagt@asu.edu, or D.
% Braghini d166353@dac.unicamp.br.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt, D. Braghini
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
% DJ - 10/19/2024: Update to use new LPI programming structure;
% DB, 24/12/2024- Update to non-coercive version

% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Pass to the 2D executive if necessary.
if PIE.dim==2
    if nargin==1
        [prog,Wo, gam] = PIETOOLS_H2_norm_2D_o(PIE);
    elseif nargin==2
        [prog,Wo, gam] = PIETOOLS_H2_norm_2D_o(PIE,settings);
    else
        [prog,Wo, gam]= PIETOOLS_H2_norm_2D_o(PIE,settings,options);
    end
    return
end
% Extract PIE operators necessary for the executive.
Top = PIE.T;        Twop = PIE.Tw;
Aop = PIE.A;        B1op = PIE.B1;
C1op = PIE.C1;

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


fprintf('\n --- Searching for H2 norm bound using the observability gramian --- \n')
% Declare an LPI program and initialize domain and opvar spaces
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
% STEP 1: declare the posopvar variable, Pop, which defines the storage
% function candidate
disp('- Declaring PI variables of the non-coercive functional...');

[prog, R1op] = poslpivar(prog, Top.dim, dd1, options1);

if override1~=1
    [prog, R2op] = poslpivar(prog, Top.dim, dd12, options12);
    Rop=R1op+R2op;
else
    Rop=R1op;
end
dim=Top.dim;
[prog,Qop] = lpivar(prog,dim,ddZ);
disp('- Declaring matrix variable...');
dimW=B1op.dim(:,2);
[prog,Wm] = poslpivar(prog,dimW);
disp('- Constructing the Equality Constraints...');
prog = lpi_eq(prog,Top'*Qop-Rop);
%prog = lpi_eq(prog,Qvar'*Top'-Rvar);
prog = lpi_eq(prog,Qop'*Top-Top'*Qop);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Using the observability gramian
disp('- Constructing the Inequality Constraints...');
ny=C1op.dim(1,1);
Dneg=[-gam*eye(ny) C1op
            C1op' Qop'*Aop+Aop'*Qop];
Dpos=[Wm B1op'*Qop
            Qop'*B1op Rop];
traceVal = trace(Wm.P);
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
prog = lpi_ineq(prog, gam-traceVal);
%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog_sol = lpisolve(prog,sos_opts); 
R = lpigetsol(prog_sol,Rop);
Q = lpigetsol(prog_sol,Qop);
W = lpigetsol(prog_sol,Wm);
if nargin<=2 || ~isfield(options,'h2')
        gam = double(lpigetsol(prog_sol,gam));
        disp('The H2 norm of the given system is upper bounded by:')
         disp(gam);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%