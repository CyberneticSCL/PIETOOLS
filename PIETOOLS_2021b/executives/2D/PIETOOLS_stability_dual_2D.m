%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_stability_dual_2D.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prog, Pop, Qop, solve_val] = PIETOOLS_stability_dual_2D(PIE,settings)
% This script executes a stability analysis for a 2D-PIE System defined
% by the PI operator representation
%   Top \dot{x}(t) = Aop x(t)
% by testing stability of the dual system
%   Top' \dot{x}(t) = Aop' x(t)
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
%   Pop:        dopvar2d object defining the LF V(x)=<x,Pop*x>. Use
%               "sosgetsol_lpivar_2d(prog,Pop)" to get the solved operator.
%   Qop:        dopvar2d object defining the right-hand side of the
%               stability LPI,
%               Qop = Top'*Pop*Aop + Aop'*Pop*Top + epneg*(Top'*Top).
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
% Initial coding DJ - 02/21/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 0: Extract LPI settings and necessary PI operators

if nargin==1
    settings = settings_PIETOOLS_light_2D;
    settings.sos_opts.simplify = 1;         % Use psimplify
    settings.eppos = [1e-4; 1e-6; 1e-6; 1e-6];    % Positivity of Lyapunov Function
    settings.epneg = 0*1e-5;                % Negativity of derivative of Lyapunov Function in both ODE and PDE state -  >0 if exponential stability desired
elseif ~isfield(settings,'is2D') || ~settings.is2D
    % Extract 2D settings.
    sos_opts = settings.sos_opts;
    settings = settings.settings_2d;
    settings.sos_opts = sos_opts;
end

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

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Extract the relevant PI operators
Aop = PIE.A;    Top = PIE.T;

% Other relevant information
vars = PIE.vars;                % Spatial variables of the PIE
dom = PIE.dom;                  % Spatial domain of the PIE
np_op = Aop.dim(:,1);           % Dimensions RxL2[x]xL2[y]xL2[x,y] of the PDE state



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Initialize the SOS program in the provided spatial variables, and
% set the Hinfty norm as an objective function value 
disp(' === Executing dual stability test === ')
prog = sosprogram(vars(:));      % Initialize the program structure
for kk=1:prod(size(vars))
    % Make sure variables are in right order (pvar stores them
    % alphabetically.
    prog.vartable(kk) = vars(kk).varname;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: declare the posopvar variable, Pop, which defines the Lyapunov 
% function candidate
disp(' - Parameterizing a positive Lyapunov operator using the specified options');

% Initialize an operator which is positive semidefinite everywhere
[prog, Pop] = poslpivar_2d(prog,np_op,dom,LF_deg,LF_opts);

% Add an additional term with psatz multiplier if requested
for j=1:length(LF_use_psatz)
    if LF_use_psatz(j)~=0
        [prog, P2op] = poslpivar_2d(prog,np_op,dom,LF_deg_psatz{j},LF_opts_psatz{j});
        Pop = Pop + P2op;
    end
end

% Ensure strict positivity of the operator
if ~all(eppos==0)
    Ip = blkdiag(eppos(1)*eye(np_op(1)),eppos(2)*eye(np_op(2)),eppos(3)*eye(np_op(3)),eppos(4)*eye(np_op(4)));
    Iop = opvar2d(Ip,[np_op,np_op],dom,vars);
    Pop = Pop + Iop;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Define the Lyapunov Inequality
%
% i.e. - Assemble the big operator
% Pheq = [ A*P*T' + T*P*A']

disp(' - Constructing the Negativity Constraint');

if epneg==0
    Qop = Top*Pop*Aop' + Aop*Pop*Top'; 
else
    % Enforce strict negativity
    Qop = Top*Pop*Aop' + Aop*Pop*Top' + epneg*(Top*Top'); 
end
% Get rid of terms that are below tolerance
ztol = 1e-12;
Qop = clean_opvar(Qop,ztol);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
%
disp(' - Enforcing the Negativity Constraint');

if use_lpi_ineq
    disp('   - Using lpi_ineq...');
    prog = lpi_ineq_2d(prog,-Qop,ineq_opts);
else
    disp('   - Using an equality constraint...');
    
    % Next, build a positive operator Qeop to enforce Qop == -Qeop;
    Qdim = Qop.dim(:,1);
    [progQ, Qeop] = poslpivar_2d(prog,Qdim,dom,eq_deg,eq_opts);
    
    toggle = 0; % Set toggle=1 to check whether the monomials in Qeop are sufficient.
    if toggle
        % Check that the parameters of Qeop indeed contain all monomials that
        % appear in the parameters of Qop
        par_indx = [2;3;4;6;7;8;11;12;16];    % Check only lower-triangular parameters
        [isgood_Qeop,isgood_Qpar,eq_deg] = checkdeg_lpi_eq_2d(Qop,Qeop,eq_deg,par_indx);

        % If Qeop is missing monomials, keep increasing the degrees until Qeop
        % has all the necessary monomials
        while ~isgood_Qeop
            warning('The specified options for the equality constraint do not allow sufficient freedom to enforce the inequality constraint. Additional monomials are being added.')

            % Construct a new positive operator Qeop with greater degrees
            [progQ, Qeop] = poslpivar_2d(prog,Qdim,dom,eq_deg,eq_opts);

            % Check that now Qeop has all the necessary monomials
            par_indx = find(~(isgood_Qpar(:)));   % Indices of parameters we still need to verify are okay
            [isgood_Qeop,isgood_Qpar,eq_deg] = checkdeg_lpi_eq_2d(Qop,Qeop,eq_deg,par_indx);
        end
    end
    prog = progQ;   % Make sure the SOS program contains the right operator Qeop

    % Introduce the psatz term.
    for j=1:length(eq_use_psatz)
        if eq_use_psatz(j)~=0
            [prog, Qe2op] = poslpivar_2d(prog,Qdim,dom, eq_deg_psatz{j},eq_opts_psatz{j});
            Qeop = Qeop+Qe2op;
        end
    end
    
    % Enforce the equality constraint Qop=-Qeop.
    prog = lpi_eq_2d(prog,Qeop+Qop);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Solve the problem, and extract the solution
%
disp(' - Solving the LPI using the specified SDP solver');
prog = sossolve(prog,sos_opts); 

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