function [prog,W, gam, R, Q] = PIETOOLS_H2_norm_2D_c_non_coercive(PIE, settings,options)
% This function solves a minimization problem to obtain the H2 norm of a linear distributed parameter
% system using the controlability gramian approach and PIEs framework. For
% the feasibility test an additional input with a assigned value to the
% norm is required.
% inputs: (mandatory)
%   PIE : PIE structure of the corresponding system
%   settings : options related to PIETOOLS
% outputs:
%   gam = computed H2 norm for SDP problem or optional input value for feasibility tests;
%   W = Observability gramian, a positive PI operator;
%  prog= sum of squares program information;
%   R = positive PI operator parameterizing the non-coercive Lyapunov 
%       functional V(v)=<v,R*v>;
%   Q = indefinite PI operator such that R=PIE.T'*Q=Q'*PIE.T;
%
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% S. Shivakumar at sshivak8@asu.edu, D. Jagt at djagt@asu.edu, or D.
% Braghini d166353@dac.unicamp.br.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  PIETOOLS Team
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
% DB, 01/17/2025: initial code

% STEP 0: Extract LPI settings and necessary PI operators
%     if nargin==1
%         settings = settings_PIETOOLS_light_2D;
%         settings.sos_opts.simplify = 1;         % Use psimplify
%         settings.eppos = [1e-4; 1e-6; 1e-6; 1e-6];    % Positivity of Lyapunov Function
%         settings.epneg = 0*1e-5;                % Negativity of derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
%     elseif ~isfield(settings,'is2D') || ~settings.is2D

% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Extract the relevant PI operators.
Top = PIE.T;        Twop = PIE.Tw;
Aop = PIE.A;        B1op = PIE.B1;
C1op = PIE.C1;

% Make sure thera are no disturbances at the boundary.
if ~(Twop==0)
    error('H2 norm LPI cannot currently be solved for systems with disturbances at the boundary');
end

% Extract 2D settings.
sos_opts = settings.sos_opts;
settings = settings.settings_2d;
settings.sos_opts = sos_opts;


% Extract options for sossolve
if ~isfield(settings,'sos_opts')
    sos_opts.simplify = 1;         % Use psimplify
    sos_opts.solver = 'sedumi';
    % % Other optional SDP solvers supported by PIETOOLS
    % settings.sos_opts.solver = 'mosek';
    % settings.sos_opts.solver = 'sdpnalplus';
    % settings.sos_opts.solver = 'sdpt3';
else
   sos_opts = settings.sos_opts; 
end

% Extract LF positivity and negativity conditions
if ~isfield(settings,'eppos')
    eppos = [1e-4;      % Positivity of Lyapunov Function with respect to real-valued states
              1e-6;1e-6;
              1e-6];    % Positivity of Lyapunov Function with respect to 2D spatially distributed states
else
    eppos = settings.eppos;
end
if ~isfield(settings,'epneg')
    epneg = 0;         % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
else
    epneg = settings.epneg;
end

% Extract settings defining operator P parameterizing LF V=<x,Px>
LF_deg = settings.LF_deg;
LF_opts = settings.LF_opts;

% Does P include a psatz term? If so, what degrees and options?
LF_use_psatz = settings.LF_use_psatz;
LF_deg_psatz = extract_psatz_deg(settings.LF_deg_psatz,LF_use_psatz);
LF_opts_psatz = extract_psatz_opts(settings.LF_opts_psatz,LF_use_psatz);

% Do we use lpi_ineq or not? Extract the appropriate options
use_lpi_ineq = settings.use_sosineq;
if use_lpi_ineq
    ineq_opts = settings.ineq_opts;
else
    eq_opts = settings.eq_opts;
    eq_deg = settings.eq_deg;
    
    eq_use_psatz = settings.eq_use_psatz;
    eq_deg_psatz = extract_psatz_deg(settings.eq_deg_psatz,eq_use_psatz);
    eq_opts_psatz = extract_psatz_opts(settings.eq_opts_psatz,eq_use_psatz);
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
    prog = lpisetobj(prog, gam);        % set gamma as objective function to minimize. This automatically set gam as a decision variable.  
else
    gam = options.h2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the storage 
% function candidate
disp('- Declaring Gramian using specified options...');

% Initialize the non-coercive operator 
[prog, Rop] = poslpivar_2d(prog, Top.dim, LF_deg, LF_opts);  
% Add an additional term with psatz multiplier if requested
for j=1:length(LF_use_psatz)
    if LF_use_psatz(j)~=0
                [prog, R2op] = poslpivar_2d(prog, Rop.dim, LF_deg_psatz{j}, LF_opts_psatz{j});
                 Rop = Rop + R2op;
    end
end

disp('- Constructing the Negativity Constraint...');

% Also declare an indefinite operator Qop so that Rop = Top*Qop.  
Qdeg = get_lpivar_degs(Rop,Top);
[prog, Qop] = lpivar_2d(prog,Top.dim,Qdeg);
prog = lpi_eq_2d(prog, Top*Qop-Rop);
    
% Finally, declare a positive operator (matrix) Wm representing the
% controllability Gramian
[prog, Wm] = poslpivar_2d(prog,C1op.dim(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STEP 2: Using the observability gramian

disp('- Constructing the Negativity Constraint...');
Iw = mat2opvar(eye(size(B1op,2)), B1op.dim(:,2), PIE.vars, PIE.dom);

Dneg = [-gam*Iw   B1op'
        B1op      Qop'*Aop'+Aop*Qop];
Dpos = [Wm           C1op*Qop
        Qop'*C1op'   Rop];

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Inequality Constraints...');

 if use_lpi_ineq
    disp('   - Using lpi_ineq...');
    prog = lpi_ineq_2d(prog,-Dneg,ineq_opts);
    prog = lpi_ineq_2d(prog,Dpos,ineq_opts);

else
    disp('   - Using an equality constraint...');
     [prog, De1op] = poslpivar_2d(prog, Dneg.dim, eq_deg, eq_opts);
     [prog, De3op] = poslpivar_2d(prog, Dpos.dim, eq_deg, eq_opts);
    
    % Introduce the psatz term.
    for j=1:length(eq_use_psatz)
        if eq_use_psatz(j)~=0
            eq_opts_psatz{j}.exclude = eq_opts_psatz{j}.exclude | eq_opts.exclude;
            eq_opts_psatz{j}.sep = eq_opts_psatz{j}.sep | eq_opts.sep;
             [prog, De2op] = poslpivar_2d(prog, Dneg.dim, eq_deg_psatz{j}, eq_opts_psatz{j});
             [prog, De4op] = poslpivar_2d(prog, Dneg.dim, eq_deg_psatz{j}, eq_opts_psatz{j});

            Deop = De1op+De2op;
            Deopp = De3op+De4op;
        end

    end
    prog = lpi_eq_2d(prog,Deop+Dneg,'symmetric'); %Dneg=-Deop
    prog = lpi_eq_2d(prog,Deopp-Dpos,'symmetric'); %Dpos=-Deopp
 end
 
% ensuring scalar inequality gam>trace
traceVal = trace(Wm.R00);
prog = lpi_ineq(prog, gam-traceVal);

% solving the lpi program
disp('- Solving the LPI using the specified SDP solver...');
prog = lpisolve(prog,sos_opts); 

W = lpigetsol(prog,Wm);
R = lpigetsol(prog,Rop);
Q = lpigetsol(prog,Qop);

if nargin<=2 || ~isfield(options,'h2')
    gam = double(lpigetsol(prog,gam));
    disp('The H2 norm of the given system is upper bounded by:')
    disp(gam);
end

end

function outcell = extract_psatz_deg(incell,use_psatz)
% Check the number of elements of incell against those of use_psatz to
% determine if a cell element has been defined for each psatz term.

if all(use_psatz==0)
    outcell = {};
    return
end
if isa(incell,'struct')
    for j=1:length(use_psatz)
        outcell{j} = incell;
    end
elseif numel(incell)==1
    for j=1:length(use_psatz)
        outcell{j} = incell{1};
    end
elseif numel(incell)==length(use_psatz)
    outcell = incell;
elseif numel(incell)>=max(use_psatz)
    outcell = incell(use_psatz);
else
    error('For each element of ''use_psatz'', a ''deg'' field should be defined.')
end

end
function outcell = extract_psatz_opts(incell,use_psatz)
% Check the number of elements of incell against those of use_psatz to
% determine if a cell element has been defined for each psatz term.

if all(use_psatz==0)
    outcell = {};
    return
end
if isa(incell,'struct')
    for j=1:length(use_psatz)
        outcell{j} = incell;
    end
elseif numel(incell)==1
    for j=1:length(use_psatz)
        outcell{j} = incell{1};
    end
elseif numel(incell)==length(use_psatz)
    outcell = incell;
elseif numel(incell)>=max(use_psatz)
    outcell = incell(use_psatz);
else
    error('For each element of ''use_psatz'', an ''opts'' field should be defined.')
end

end