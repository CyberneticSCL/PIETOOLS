function [prog,Wc, gam] = PIETOOLS_H2_norm_2D_c(PIE, settings,options)
% This function solves a minimization problem to obtain the H2 norm of a linear distributed parameter
% system using the controllability gramian approach and PIEs framework. For
% the feasibility test an additional input with a assigned value to the
% norm is required.
% inputs: (mandatory)
%   PIE : PIE structure of the corresponding system
%   settings : options related to PIETOOLS
% outputs:
%   gam = computed H2 norm for SDP problem or optional input value for feasibility tests;
%   Wc= Controllability gramian, a positive PI operator;
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
% DJ - 10/20/2024: Update to use new LPI programming structure;
% DJ, 01/17/2025: Bugfix Pop --> Wop;

% STEP 0: Extract LPI settings and necessary PI operators


% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Extract the relevant PI operators.
Top = PIE.T;        Twop = PIE.Tw;
Aop = PIE.A;        Bwop = PIE.B1;
Czop = PIE.C1;

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
    % prog = lpidecvar(prog, gam);        % set gam = gamma as decision variable
    % prog = lpi_ineq(prog, gam);         % enforce gamma>=0
    prog = lpisetobj(prog, gam);        % set gamma as objective function to minimize
else
    gam = options.h2;
end
%
% Alternatively, the above 3 commands may be commented and a specific gain
% test specified by defining a specific desired value of gamma. This
% results in a feasibility test instead of an optimization problem.
% gamma = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: declare the posopvar variable, Pop, which defines the storage 
% function candidate
disp('- Declaring Gramian using specified options...');

% Initialize an operator which is positive semidefinite everywhere
[prog, Wop] = poslpivar_2d(prog, Top.dim, LF_deg, LF_opts);                 % DJ, 01/17/2025: Pop --> Wop;
%[prog, P1op] = poslpivar(prog, [nx1 ,nx2],X,dd1,options1);

% Add an additional term with psatz multiplier if requested
for j=1:length(LF_use_psatz)
    if LF_use_psatz(j)~=0
        [prog, P2op] = poslpivar_2d(prog,Wop.dim, LF_deg_psatz{j}, LF_opts_psatz{j});
        Wop = Wop + P2op;
    end
end

% Ensure strict positivity of the operator
if ~all(eppos==0)
    np_op = Wop.dim(:,1);           % Dimensions RxL2[x]xL2[y]xL2[x,y] of the PDE state
    Ip = blkdiag(eppos(1)*eye(np_op(1)),zeros(np_op(2)),zeros(np_op(3)),eppos(4)*eye(np_op(4)));
    Iop = opvar2d(Ip, Wop.dim, PIE.dom, PIE.vars);
    Wop = Wop + Iop;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Using the controlability gramian

disp('- Constructing the Negativity Constraint...');

Dop =  (Aop*Wop)*Top'+Top*(Wop*Aop')+Bwop*Bwop';
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp('- Enforcing the Inequality Constraints...');

 if use_lpi_ineq
    disp('   - Using lpi_ineq...');
    prog = lpi_ineq_2d(prog,-Dop,ineq_opts);
else
    disp('   - Using an equality constraint...');
    [prog, Deop] = poslpivar_2d(prog, Dop.dim, eq_deg, eq_opts);
    
    % Introduce the psatz term.
    for j=1:length(eq_use_psatz)
        if eq_use_psatz(j)~=0
            eq_opts_psatz{j}.exclude = eq_opts_psatz{j}.exclude | eq_opts.exclude;
            eq_opts_psatz{j}.sep = eq_opts_psatz{j}.sep | eq_opts.sep;
            [prog, De2op] = poslpivar_2d(prog, Dop.dim, eq_deg_psatz{j}, eq_opts_psatz{j});
            Deop = Deop+De2op;
        end
    end
    prog = lpi_eq_2d(prog,Deop+Dop,'symmetric'); %Dop=-Deop
 end

tempObj = Czop*Wop*Czop';
tempMat = tempObj.R00;
traceVal=0;
for idx = 1:size(tempMat,1)
    traceVal = traceVal+tempMat(idx,idx);
end
% traceVal>=gam
prog = lpi_ineq(prog, gam-traceVal);


%solving the sos program
disp('- Solving the LPI using the specified SDP solver...');
prog = lpisolve(prog,sos_opts); 

Wc = lpigetsol(prog,Wop);

if nargin<=2 || ~isfield(options,'h2')
    gam = sqrt(double(lpigetsol(prog,gam)));
    disp('The H2 norm of the given system is upper bounded by:')
    disp(gam);
end
%% set outputs
%     output.h2=gam;
%     output.W=P;
%     output.prog=prog;

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