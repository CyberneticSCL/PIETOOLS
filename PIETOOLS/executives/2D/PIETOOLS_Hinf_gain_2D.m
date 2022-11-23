%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_Hinf_gain_2D.m     PIETOOLS 2022a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prog, Pop, gam, solve_val] = PIETOOLS_Hinf_gain_2D(PIE, settings, gain)
% This function use the KYP lemma in primal form to compute an upper bound
% on the H-infty gain a 2D-PIE system of the form
%
% T_op \dot{x}(t) = A_op  x(t) + Bw_op  w(t)
%            z(t) = Cz_op x(t) + Dzw_op w(t)
%
% If any other parts of the PIE are present, these are ignored. T_op, A_op,
% Bw_op, Cz_op, and Dzw_op must be properly defined for the script to function.
%
% INPUTS:
%   PIE:        A structure defining the PIE for which to estimate the
%               H_infty gain. Must contain at least opvar2d objects T, Tw,
%               A, Bw, Cz, Dzw, as well as a 2x2 array "dom" defining the
%               domain, and a 2x2 pvar object "vars" describing the spatial
%               variables. All other fields will be ignored.
%   settings:   A structure specifying accuracy/complexity conditions of
%               the LPI; see the "settings" files for a full list of
%               included settings. 
%   gain:       (optional) An upper bound on the H_infty gain to verify
%               If no gain is provided, the solver will look for a smallest
%               upper bound on this gain.
%
% OUTPUTS:
%   prog:       SOS program structure associated to the (solved) LPI
%   Pop:        dopvar2d object defining the LF V(x)=<x,Pop*x>. Use
%               "sosgetsol_lpivar_2d(prog,Pop)" to get the solved operator.
%   gam:        Upper bound on the H_infty gain found by the solver.
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
elseif nargin==2 
    gain = 0;
end
if ~isfield(settings,'is2D') || ~settings.is2D
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
    % settings.sos_opts.solver ='mosek';
    % settings.sos_opts.solver='sdpnalplus';
    % settings.sos_opts.solver='sdpt3';
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

% Do we use bisection to estimate the Hinf norm?
if isfield(settings,'use_bisect')
    use_bisect = settings.use_bisect;
    if use_bisect
        if isfield(settings.bisect_opts,'Nruns')
            Nruns = settings.bisect_opts.Nruns; % Number of times to bisect
        else
            Nruns = 10;
        end
        if isfield(settings.bisect_opts,'min')
            gam_min = settings.bisect_opts.min; % Smallest gain to test
        else
            gam_min = 0;
        end
        if isfield(settings.bisect_opts,'max')
            gam_max = settings.bisect_opts.max; % Largest gain to test
        else
            gam_max = 1000;
        end
        if isfield(settings.bisect_opts,'start') && gain==0
            gain = settings.bisect_opts.start;  % Starting gain
        else
            gain = 0.5*gam_min + 0.5*gam_max;
        end
    end
else
    use_bisect = false;
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Extract the relevant PI operators
Twop = PIE.Tw;    
Aop = PIE.A;
Top = PIE.T;
Bop = PIE.B1;
Cop = PIE.C1;
Dop = PIE.D11;
if ~(Twop==0) && ~(Bop==0)
    error('The PIE takes both the input w and its derivative dw/dt; LPI based H_infty gain analysis is currently not supported')
end

% Other relevant information
vars = PIE.vars;                % Spatial variables of the PIE
dom = PIE.dom;                  % Spatial domain of the PIE
np_op = Aop.dim(:,1);           % Dimensions RxL2[x]xL2[y]xL2[x,y] of the PDE state
nw_op = Bop.dim(:,2);           % Dimensions of the input w
nz_op = Cop.dim(:,1);           % Dimensions of the output z
nw = sum(nw_op);                % Total number of disturbance signal
nz = sum(nz_op);                % Total number of regulated output signals



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Initialize the SOS program in the provided spatial variables, and
% set the Hinfty norm as an objective function value 

prog = sosprogram(vars(:));         % Initialize the program structure
for kk=1:prod(size(vars))
    % Make sure variables are in right order (pvar stores them
    % alphabetically.
    prog.vartable(kk) = vars(kk).varname;
end
if gain==0
    fprintf('\n --- Searching for an Hinf gain bound using the primal KYP lemma --- \n')
    % If no gain is provided, include the gain as a decision variable, and
    % look for a minimum
    dpvar gam;
    prog = sosdecvar(prog, gam);    % This sets gam = gamma as decision var
    prog = sossetobj(prog, gam);    % This minimizes gamma, comment for feasibility test
else
    fprintf(['\n --- Testing the Hinf gain bound gam = ',num2str(gain),' using the primal KYP lemma --- \n'])
    % If a gain is provided, just check the LPI is feasible with this gain
    gam = gain;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Declare the posopvar variable, Pop, which defines the storage 
% function candidate
disp('- Declaring a positive Lyapunov operator variable using the specified options...');

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
% STEP 3: Define the KYP operator
disp('- Constructing the negativity constraint...');

% Define identity PI operators
Iwop = opvar2d(eye(nw),[nw_op,nw_op],dom,vars);
Izop = opvar2d(eye(nz),[nz_op,nz_op],dom,vars);

% Assemble the KYP operator
if epneg==0
    Qop = [-gam*Iwop+Twop'*(Pop*Bop)+(Pop*Bop)'*Twop,   Dop',        (Pop*Bop)'*Top+Twop'*(Pop*Aop);
           Dop,                                         -gam*Izop,   Cop;
           Top'*(Pop*Bop)+(Pop*Aop)'*Twop,              Cop',        (Pop*Aop)'*Top+Top'*(Pop*Aop)];
else
    % Ensure strict negativity
    Qop = [-gam*Iwop+Twop'*(Pop*Bop)+(Pop*Bop)'*Twop,   Dop',        (Pop*Bop)'*Top+Twop'*(Pop*Aop);
           Dop,                                         -gam*Izop,   Cop;
           Top'*(Pop*Bop)+(Pop*Aop)'*Twop,              Cop',        (Pop*Aop)'*Top+Top'*(Pop*Aop) + epneg*(Top'*Top)];
end
% Get rid of terms that are below tolerance
ztol = 1e-12;
Qop = clean_opvar(Qop,ztol);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Impose Negativity Constraint. There are two methods, depending on
% the options chosen
disp('- Enforcing the negativity constraint...');

if use_lpi_ineq
    disp('  - Using lpi_ineq...');
    prog = lpi_ineq_2d(prog,-Qop,ineq_opts);
else
    disp('  - Using an equality constraint...');
    
    % First, check if we can exclude any monomials
    tol = ztol*1e-1;
    eq_opts_psatz_exclude = zeros(1,16);
    if (isa(Qop.R22{1,1},'double') && max(max(Qop.R22{1,1}))<tol) || all(max(max(Qop.R22{1,1}.C))<tol)
        eq_opts.exclude(8) = 1;
        eq_opts_psatz_exclude(8) = 1;
    end
    if (isa(Qop.R22{2,1},'double') && max(max(Qop.R22{2,1}))<tol) || all(max(max(Qop.R22{2,1}.C))<tol)
        eq_opts.exclude(9) = 1;
        eq_opts_psatz_exclude(9) = 1;
    end
    if (isa(Qop.R22{3,1},'double') && max(max(Qop.R22{3,1}))<tol) || all(max(max(Qop.R22{3,1}.C))<tol)
        eq_opts.exclude(10) = 1;
        eq_opts_psatz_exclude(10) = 1;
    end
    if (isa(Qop.R22{1,2},'double') && max(max(Qop.R22{1,2}))<tol) || all(max(max(Qop.R22{1,2}.C))<tol)
        eq_opts.exclude(11) = 1;
        eq_opts_psatz_exclude(11) = 1;
    end
    if (isa(Qop.R22{1,3},'double') && max(max(Qop.R22{1,3}))<tol) || all(max(max(Qop.R22{1,3}.C))<tol)
        eq_opts.exclude(12) = 1;
        eq_opts_psatz_exclude(12) = 1;
    end
    
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
    prog = progQ;    % Make sure the SOS program contains the right operator Qeop
    
    % Introduce the psatz term.
    for j=1:length(eq_use_psatz)
        if eq_use_psatz(j)~=0
            eq_opts_psatz{j}.exclude(8:12) = eq_opts_psatz_exclude(8:12);
            [prog, Qe2op] = poslpivar_2d(prog,Qdim,dom, eq_deg_psatz{j},eq_opts_psatz{j});
            Qeop = Qeop+Qe2op;
        end
    end
    
    % Enforce the equality constraint.
    prog = lpi_eq_2d(prog,Qeop+Qop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Solve the problem, and check the solution
disp('- Solving the LPI using the specified SDP solver...');

% Solve the sos program
prog = sossolve(prog,sos_opts); 

% Check the results
if norm(prog.solinfo.info.feasratio-1)<=.3 && ~prog.solinfo.info.numerr && ~prog.solinfo.info.pinf && ~prog.solinfo.info.dinf
    if gain==0
        gam = double(sosgetsol(prog,gam));
    end
    disp('The H-infty norm of the given system is upper bounded by:')
    disp(gam); % check the Hinf norm, if the problem was solved successfully
    solve_val = 1;
elseif norm(prog.solinfo.info.feasratio-1)<=.3 && prog.solinfo.info.numerr && ~prog.solinfo.info.pinf && ~prog.solinfo.info.dinf
    if gain==0
        gam = double(sosgetsol(prog,gam));
    end
    disp('The system of equations was successfully solved. However, Double-check the precision.')
    disp('The H-infty norm of the given system is upper bounded by:')
    disp(gam);
    solve_val = 0.5;
elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || norm(prog.solinfo.info.feasratio+1)<=.1
    disp('The system of equations was not solved.')
    if gain==0
        gam = double(sosgetsol(prog,gam));
    end
    disp(gam);
    %gam = Inf;
    solve_val = 0;
else
    disp('Unable to definitively determine feasibility.')
    if gain==0
        gam = double(sosgetsol(prog,gam));
    end
    disp(gam);
    %gam = Inf;
    solve_val = 0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: Repeat?
% If we're bisecting, continue with a new value of gamma
if use_bisect
    % Store the results as arrays
    gam_arr = nan*ones(1,Nruns);
    gam_arr(1) = gam;
    solve_arr = nan*ones(1,Nruns);
    solve_arr(1) = solve_val;
    gam_old = gain;
    
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    for k=2:Nruns
        % Update the value of gam
        if solve_val==1
            gam_max = gam_old;
            gam = 0.5*(gam_min+gam_max);
        elseif solve_val==0.5
            gam_max = gam_old;
            gam = 0.25*gam_min + 0.75*gam_max;
        elseif solve_val==0
            gam_min = gam_old;
            gam = 0.5*(gam_min+gam_max);
        end
        
        fprintf(['\n -- Testing the Hinf gain bound gam = ',num2str(gam),' -- \n'])
        
        % STEP 3b: Define the KYP operator
        disp('- Updating the negativity constraint...');
        gam_diff = gam - gam_old;
        Oop = opvar2d([],[np_op,np_op],dom,vars);   % Zero operator
        Qop_diff = blkdiag(-gam_diff*Iwop,-gam_diff*Izop,Oop);
        Qop = Qop + Qop_diff;
        
        % STEP 4b: Impose Negativity Constraint. There are two methods,
        %  depending on the options chosen
        disp('- Enforcing the new negativity constraint...');
        
        if use_lpi_ineq
            disp('  - Using lpi_ineq...');
            prog = lpi_ineq_2d(prog,-Qop,ineq_opts);
        else
            disp('  - Using an equality constraint...');
            [prog, Qeop] = poslpivar_2d(prog,np_op+nw_op+nz_op,dom,eq_deg,eq_opts);
            for j=1:length(eq_use_psatz)
                if eq_use_psatz(j)~=0
                    [prog, Qe2op] = poslpivar_2d(prog,np_op+nw_op+nz_op,dom, eq_deg_psatz{j},eq_opts_psatz{j});
                    Qeop = Qeop+Qe2op;
                end
            end
            prog = lpi_eq_2d(prog,Qeop+Qop); %Dop=-Deop
        end
        
        % STEP 5b: Solve the problem, and check the solution
        disp('- Solving the LPI using the specified SDP solver...');
        
        % Solve the sos program
        prog = sossolve(prog,sos_opts);
        
        % Check the results
        if norm(prog.solinfo.info.feasratio-1)<=.3 && ~prog.solinfo.info.numerr && ~prog.solinfo.info.pinf && ~prog.solinfo.info.dinf
            disp('The H-infty norm of the given system is upper bounded by:')
            disp(gam); % check the Hinf norm, if the problem was solved successfully
            solve_val = 1;
        elseif norm(prog.solinfo.info.feasratio-1)<=.3 && prog.solinfo.info.numerr && ~prog.solinfo.info.pinf && ~prog.solinfo.info.dinf
            disp('The system of equations was successfully solved. However, Double-check the precision.')
            disp('The H-infty norm of the given system is upper bounded by:')
            disp(gam);
            solve_val = 0.5;
        elseif prog.solinfo.info.pinf || prog.solinfo.info.dinf || norm(prog.solinfo.info.feasratio+1)<=.1
            disp('The system of equations was not solved.')
            solve_val = 0;
        else
            disp('Unable to definitively determine feasibility.')
            solve_val = 0;
        end
        
        % Store the results
        gam_arr(k) = gam;
        solve_arr(k) = solve_val;
        gam_old = gam;
        
    end
    
    % Return the full arrays as output
    gam = gam_arr;
    solve_val = solve_arr;
    
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