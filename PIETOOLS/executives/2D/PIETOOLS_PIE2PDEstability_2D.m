%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_PIE2PDEstability_2D.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prog, Qop, Dop, solve_val] = PIETOOLS_PIE2PDEstability_2D(PIE,settings)
% This script executes a PIE to PDE stability analysis for a 2D-PIE System
% defined by the PI operator representation
%       Top \dot{x}(t)=Aop x(t)
% testing whether solutions satisfy
%       ||Top*x(t)||_{L2} <= M*e^{-epneg t} ||x(0)||_{L2}
%
% If any other parts of the PIE are present, these are ignored. Both Top
% and Aop must be properly defined for the script to function.
%
% INPUTS:
%   PIE:        A structure defining the PIE for which to test stability. 
%               Must contain at least opvar2d objects Top and Aop, as well 
%               as a 2x2 array "dom" defining the domain, and a 2x2 pvar 
%               object "vars" describing the spatial variables. All other 
%               fields will be ignored.
%   settings:   A structure specifying accuracy/complexity conditions of
%               the LPI; see the "settings" files for a full list of
%               included settings. 
%
% OUTPUTS:
%   prog:       SOS program structure associated to the (solved) LPI.
%   Qop:        dopvar2d object defining the LF V(x)=<x,Pop*x>. Use
%               "sosgetsol_lpivar_2d(prog,Pop)" to get the solved operator.
%   Dop:        dopvar2d object defining the right-hand side of the
%               stability LPI,
%               Dop = Qop'*Aop + Aop'*Qop + epneg*(Top'*Qop).
%   solve_val:  Value specifying whether the problem is feasible, based on
%               outputs of SOS program. Value is 1 if the a feasible
%               solution is found, 0.5 if it a feasible solution is only
%               found up to reduced accuracy, and 0 if no feasible solution 
%               is found.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ - 04/05/2026: Initial coding;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0: Extract LPI settings and necessary PI operators

% Check if the PIE is properly specified.
if ~isa(PIE,'pie_struct')
    error('The PIE for which to run the executive should be specified as object of type ''pie_struct''.')
else
    PIE = initialize(PIE);
end
% Extract the relevant PI operators.
Aop = PIE.A;    Top = PIE.T;

if nargin==1
    settings = settings_PIETOOLS_light_2D;
    settings.sos_opts.simplify = 1;         % Use psimplify
    settings.eppos = 1e-2*ones(4,1);    % Positivity of Lyapunov Function
    settings.epneg = 0*1e-5;                % Negativity of derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
elseif ~isfield(settings,'is2D') || ~settings.is2D
    % Extract 2D settings.
    sos_opts = settings.sos_opts;
    settings = settings.settings_2d;
    settings.sos_opts = sos_opts;
end

% Set tolerance below which coefficients are assumed to be zero
ztol = 1e-12;

% Extract options for sossolve
if ~isfield(settings,'sos_opts')
    sos_opts.simplify = 1;         % Use psimplify
    %sos_opts.solver = 'sedumi';
    % % Other optional SDP solvers supported by PIETOOLS
    % settings.sos_opts.solver = 'mosek';
    % settings.sos_opts.solver = 'sdpnalplus';
    % settings.sos_opts.solver = 'sdpt3';
else
   sos_opts = settings.sos_opts; 
end

% Extract LF positivity and negativity conditions
if ~isfield(settings,'eppos')
    eppos = 1e-2*ones(4,1);      % Positivity of Lyapunov Function
else
    eppos = settings.eppos;
    if isscalar(eppos)
        eppos = eppos*ones(4,1);
    end
end
if ~isfield(settings,'epneg')
    epneg = 0;         % Negativity of Derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
else
    epneg = settings.epneg;
end

% Extract settings defining operator Q parameterizing LF V=<x,Rx>
Qop_deg = settings.LF_deg;
Qop_opts = rmfield(settings.LF_opts,'psatz');
Rop_deg = Qop_deg;

% Do we enforce R>=0 using lpi_ineq, or using manual construction?
Rop_use_ineq = settings.use_sosineq;

% Does R include a psatz term? If so, what degrees and options?
Rop_use_psatz = settings.LF_use_psatz;
Rop_deg_psatz = extract_psatz_deg(settings.LF_deg_psatz,Rop_use_psatz);
Rop_opts_psatz = extract_psatz_opts(settings.LF_opts_psatz,Rop_use_psatz);

% Do we use lpi_ineq or not? Extract the appropriate options
use_lpi_ineq = settings.use_sosineq;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Initialize the SOS program in the provided spatial variables, and
% set the Hinfty norm as an objective function value 
disp(' === Executing primal stability test === ')
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);      % Initialize the program structure



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: declare the posopvar variable, Rop, which defines the Lyapunov 
% function candidate
disp(' - Parameterizing a coercive Lyapunov functional with respect to the PDE state');


% % First, declare an indefinite operator Qop, representing Qop=Pop*Top.
[prog, Qop] = lpivar_2d(prog, Top.dim, Qop_deg, Qop_opts);
TQop = Top'*Qop;

% % Then, enforce Top'*Qop >= 0.
if Rop_use_ineq
    % Use 'lpi_ineq' to generate Rop
    Rop_opts.psatz = Rop_use_psatz;
    [prog, ~] = lpi_ineq_2d(prog, TQop, Rop_opts);
else
    % Manually generate Rop
    Rop_opts.psatz = 0;
    Rop_opts.exclude = zeros(1,16);
    Rop_opts.sep = zeros(6,1);
    Rop_opts = get_eq_opts_2D(TQop,Rop_opts,ztol);
    [prog, Rop] = poslpivar_2d(prog,Top.dim, Rop_deg, Rop_opts);

    % Introduce the psatz term.
    for j=1:length(Rop_use_psatz)
        if Rop_use_psatz(j)~=0
            Rop_deg_j = Rop_deg_psatz{j};
            Rop_opts_j = Rop_opts_psatz{j};
            Rop_opts_j.psatz = Rop_use_psatz(j);
            [prog, R2op] = poslpivar_2d(prog, Rop.dim, Rop_deg_j, Rop_opts_j);
            Rop = Rop+R2op;
        end
    end

    % Enforce the equality constraint Top'*Qop = Rop.
    prog = lpi_eq_2d(prog, TQop-Rop);
end

% % Finally, ensure Top'*Qop >= eppos*Top'*Top by setting
% %     Qop --> Qop + eppos*Top
np_op = Top.dim(:,1);           % Dimensions RxL2[x]xL2[y]xL2[x,y] of the PDE state
Ip = blkdiag(eppos(1)*eye(np_op(1)),eppos(2)*eye(np_op(2)),eppos(3)*eye(np_op(3)),eppos(4)*eye(np_op(4)));
Iop = mat2opvar(Ip,Top.dim,[Top.var1,Top.var2],Top.I);
Qop = Qop + Iop*Top;
TQop = Top'*Qop;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Compute the Lie derivative of the Lyapunov functional along the
% PIE
%
%   d/dt V(v(t)) = d/dt < v(t), Qop'*Top*v(t)> 
%                = < v(t), [Qop'*Aop + Aop'*Qop]*v(t)>

disp(' - Computing the Lie derivative of the Lyapunov functional along the PIE');

AQop = Aop'*Qop;
Dop = AQop' + AQop; 
if epneg~=0
    % Enforce strict negativity
    Dop = Dop + 2*epneg*TQop; 
end
% Get rid of terms that are below tolerance
Dop = clean_opvar(Dop,ztol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp(' - Enforcing the Negativity Constraint');

if use_lpi_ineq
    disp('   - Using lpi_ineq...');
    prog = lpi_ineq_2d(prog,-Dop,settings.ineq_opts);
else
    disp('   - Using an equality constraint...');
    
    % % Next, build a positive operator Qeop to enforce Qop == -Qeop;
    eq_deg = settings.eq_deg;
    eq_opts.psatz = 0;
    eq_opts.sep = zeros(6,1);
    eq_opts.exclude = zeros(1,16);
    eq_opts = get_eq_opts_2D(Dop,eq_opts,ztol);
    [prog, Deop] = poslpivar_2d(prog,Dop.dim,eq_deg,eq_opts);
    
    % Introduce the psatz term.
    diff_use_psatz = settings.eq_use_psatz;
    for j=1:length(diff_use_psatz)
        if diff_use_psatz(j)~=0
            eq_deg_j = settings.eq_deg_psatz{j};
            eq_opts_j = eq_opts;
            eq_opts_j.psatz = diff_use_psatz(j);
            [prog, De2op] = poslpivar_2d(prog, Dop.dim, eq_deg_j, eq_opts_j);
            Deop = Deop+De2op;
        end
    end
    
    % Enforce the equality constraint Qop=-Qeop.
    prog = lpi_eq_2d(prog,Deop+Dop,'symmetric');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Solve the problem, and extract the solution
%
disp(' - Solving the LPI using the specified SDP solver');
prog = lpisolve(prog,sos_opts); 

% Conclusion:
if norm(prog.solinfo.info.feasratio-1)<=.3 && ~prog.solinfo.info.numerr && ~prog.solinfo.info.pinf && ~prog.solinfo.info.dinf
    disp('The System of equations was successfully solved.')
    solve_val = 1;
elseif norm(prog.solinfo.info.feasratio-1)<=.3 && prog.solinfo.info.numerr && ~prog.solinfo.info.pinf && ~prog.solinfo.info.dinf
    disp('The System of equations was successfully solved. However, Double-check the precision.')
    solve_val = 0.5;
elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || norm(prog.solinfo.info.feasratio+1)<=.1
    disp('The System of equations was not solved.')
    solve_val = 0;
else
    disp('Unable to definitively determine feasibility. Numerical errors dominating or at the limit of stability.')
    solve_val = 0;
end

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
elseif isscalar(incell)
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


%%
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
elseif isscalar(incell)
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