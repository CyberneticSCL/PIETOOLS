%
% Test the SOS stability test on the 2D linear PDEs in the PDE library
%
% This script loops over the 2D stability examples (33-44) of the library
% `examples_PDE_library_PIETOOLS_nouinput`, and for each one repeats the
% steps of `Test_Linear_Examples`: it converts the PDE to a PIE, declares a
% Lyapunov functional V with `piesos_LF`, and enforces dV <= 0 with
% `piesos_ineq`. Each program is solved, and its status is compared to the
% expected stability of the example. With u=T*x the PDE state, V has one
% of two forms, set by the field lf of the settings:
%   'pde':  V(x) = <Zop*Z(u),Pmat*Zop*Z(u)> + eppos*||u||^2, Pmat>=0, the
%           form of the stock executive `PIETOOLS_stability_2D`;
%   'pie':  V(x) = <Zop*Z(x),Pmat*Zop*Z(u)>, Pmat symmetric, with
%           V >= eppos*||u||^2 enforced here, like the LPI of Cor. 35 of
%           Jagt & Peet (2025) at rate 0. The energy of wave-type systems
%           involves derivatives of u, so is of this form, but is not
%           V(u) with operators on u, which probably is why 'pde' fails
%           for 36 and 38 (inferred, as for 16 and 18 of the 1D tests).
%
% The cost of the SDP explodes with the degrees in two spatial variables,
% so the defaults are the smallest with which the 2D reaction-diffusion
% equation (34) is certified up to 0.85 of its stability limit: Zop of
% degree 1 (def.pdegs) and monomials of degree 2 in `piesos_ineq`
% (def.ineq_degs, psatz 1). With the default degrees of `piesos_ineq`, a
% single SDP takes more than 15 minutes, and with pdegs=2 more than 5.5 GB.
%
% Examples included (linear, stability, certifiable at these degrees):
%   33  transport, c1,c2<0
%   34  reaction-diffusion, stable for lam<=2*pi^2
%   35  advection-diffusion, stable for all b1,b2
%   37  parabolic, Dirichlet-Neumann, stable for a<=pi^2/2
%   36  telegraph equation, and 38 wave equation, with lf='pie' (with
%       lf='pde' both are infeasible, see STATUS). 38 is undamped, so
%       only marginally stable (dV=0). 36 is also run with negative
%       damping, lam=-1, which is unstable, as a check on the 'pie' form
%   40  heat equation coupled to an ODE, stable for K<-1
%   42, 43  heat equation with delay, known failures (see STATUS)
% The parameters lam of 34 and a of 37 are set to fractions of their limit,
% in column 5 of ex_tab. Below the limit the program is expected feasible,
% beyond it infeasible (a certificate for an unstable system would be a
% bug). Closer to the limit the stable system may be infeasible at these
% degrees, so the outcome is not prescribed; the report gives the largest
% feasible fraction. Example 40 is also run beyond the stability limit of
% the ODE (K=-0.5). The default of 43 has +a1*x(t-tau) in the code, but
% -a1*x(t-tau) in the description of the library (stable for tau<1.0347):
% a1=-1 is stable by that description, and is run as well.
%
% Examples not included (see STATUS for the measurements):
% - 39 (Navier-Stokes): 3 states, not run (36 and 38 with 2 states take
%   4 to 5 GB with lf='pie');
% - 41, 44 (heat equation with delayed boundary input, wave equation with
%   boundary delay): infeasible, also with the classical test (`convert`
%   and `lpiscript(PIE,'stability','heavy')`);
% - 45-52: Hinf-gain, observer, H2 and controller problems.
%
% The user may change the list ex_list of examples to run (set it before
% running the script), the default settings `def` (also by a struct
% def_override), or a row of table `ex_tab`, with columns
%   {number, expected, settings, arguments, fraction},
% where expected is true (stable, must be feasible), false (unstable,
% must be infeasible) or NaN (not prescribed), settings is a struct
% with any of the fields of `def`, overriding the defaults, arguments a
% cell of extra inputs for the library, e.g. {'lam=10;'}, and fraction the
% fraction of the stability limit at which the parameter is set (NaN if
% not a sweep). For lf='pie', set ineq_degs to struct('in',d,'out',p), with
% d the degree of the dummy variables and p of the primary variables of
% the operators in `piesos_ineq`.
%
% STATUS (10/22/2026, MOSEK 11.2.2, default settings, 22 rows, 30 min in
% total): 19 rows give the expected result, 3 rows fail with an error.
% Sizes m x n are of the SDP (rows of y, x); times are for one example,
% peaks are of a fresh Matlab for one example (the report gives the peak
% of the whole process so far); on battery, times may be up to 2x too
% long:
% - 33 transport: feasible for c1=c2=-1 and for c1=-2, c2=-0.5 (20 s,
%   2282x92944, 1.7 GB);
% - 34 reaction-diffusion (lam limit 2*pi^2=19.74): feasible at 0.5, 0.75
%   and 0.85 of the limit (feasratio 0.82, 0.90, 0.74) for psatz=1, 
%   infeasible (pinf=1) at 0.9, 0.95, 0.99 and 1.05. Works for 0.95 with
%   psatz=2;
% - 35 advection-diffusion: feasible for the library parameters and for
%   b1=2, b2=4 (70-85 s, 4688x92944, 3.4 GB);
% - 37 parabolic (a limit pi^2/2=4.93): feasible at 0.5, 0.75 and 0.99
%   of the limit (also 0.9, 0.95), infeasible at 1.05 (and 1.1). Largest
%   feasible fraction 0.99 (40 s, 4002x92944, 2.2 GB);
% - 36 telegraph and 38 wave, with lf='pie' and (d,p)=(1,2) in piesos_ineq:
%   36 feasible (feasratio 1.11, 170-270 s, 35260x265536, 4.9 GB), and
%   infeasible (pinf=1) with lam=-1; 38 feasible (0.875, 145-170 s, 3.6
%   GB; undamped, so dV=0), and infeasible with c1=c2=-1 (measured outside
%   the table). With lf='pde' both are infeasible (pinf=1; ineq degs 2: 36
%   650 s, 17696x371776, 7.5 GB; 38 865 s, 6.3 GB). The degrees (1,2) were
%   chosen on a telegraph variant with Dirichlet-Neumann conditions (not in
%   the library): ineq degs 1 and (2,1) are infeasible (pinf=1), (1,2) is
%   feasible (5.0, 170 s, 4.3 GB), and (2,2) gives 986k x 36.6k, killed at
%   7 GB (also for 36);
% - 40 heat equation with ODE: feasible for K=-2 and K=-1.25 (feasratio
%   1.9, 1.3), infeasible for K=-0.5, with the ODE unstable (75 s,
%   4108x98432, 3.0 GB). The library sets K=k*eye before evaluating its
%   arguments, so 'k=..;' has no effect, use 'K=..;'.
% Known failures:
% - 42, 43 (a1=1 and a1=-1): error "Input function space of the functional
%   must match output space of the operator" in innerprod(u,u) with
%   u=PIE.T*x (stack: mtimes_tensopvar, mtimes, int, innerprod), after
%   piesos_LF has run. States depend on s1 and on (s1,s2); innerprod(x,x)
%   works, and so does innerprod(u,u) for 44, with states in s1 and s2.
% Standard LPI programming for 36 and 38 (probes, not in the table): the
% stock 2D executive (`lpiscript` 'stability', 'stripped', with V=<Tx,P Tx>
% as lf='pde') is infeasible for both (pinf=1; 36: 530 s, 5.1 GB; 38: 830
% s, 5.1 GB), but it is also infeasible for 34 at lam=9.87 (numerr=2), as is
% 'light', and 'light' with SeDuMi gives numerr=1 there and pinf=1 at 14.8,
% so it is a weak reference. The LPI of Cor. 35 at rate 0 (`lpivar` default
% degrees, P'*T=T'*P, T'*P-ep^2*T'*T>=0, -(P'*A+A'*P)>=0, ep=0.1, see
% jp_wave_plate_rates) is infeasible for 36 and 38 with the Dirichlet
% conditions of the library (36 also with psatz 3-6, `lpivar` degree 2 and
% ep=0.01; 20-160 s, up to 4.7 GB), but feasible for the telegraph with
% Dirichlet-Neumann conditions (1.31, 26 s) and for the damped wave of the
% paper (kap=1; 1.33, 50 s), so the probe is set up correctly. Why this
% LPI fails and lf='pie' succeeds for 36 and 38 was not investigated.
% Not included, measured with the same settings:
% - 39: not run, as it has 3 states;
% - 41, 44: infeasible also for pdegs=2 and ineq degs=3, and for the
%   classical test ('heavy'; 41: feasratio 0.495, numerr=2; 44: pinf=1).
%   For 41, `convert_PIETOOLS_PDE_nonlinear` merges s1 with the delay
%   variable, so the PIE has one spatial variable.
% By hand, the default of 43 (+a1*x(t-tau)) has mode 1 with characteristic
% equation s=0.9+exp(-s*tau), with a positive root, so it is unstable,
% unlike the library description (not verified numerically).
%
% DJ (CC), 10/22/2026: Initial coding;

%%%% 0. Settings

% Defaults for the Lyapunov functional and the inequality, see
% Test_Linear_Examples
def.lf = 'pde';         % form of V: 'pde' for <Z(u),P Z(u)>, 'pie' for <Z(x),P Z(u)>
def.d = 1;              % degree of the monomial basis Z(u), doubled in V
def.pdegs = 1;          % maximal monomial degree of Zop
def.exclude = [0;0;0];  % operators of Zop to exclude (mult, lower, upper)
def.eppos = 1e-4;       % coefficient of the positivity term eppos*||u||^2
def.psatz = 2;          % order of the multipliers in piesos_ineq
def.ineq_degs = 2;      % degrees of the monomials in piesos_ineq (default if [])
def.solver = '';        % SDP solver passed to piesos_solve ('' for the default)
def.feas_tol = 0.5;     % minimal feasratio to accept a solution
if exist('def_override','var')
    % Optionally override the defaults by a struct def_override
    for fn = fieldnames(def_override)'
        def.(fn{1}) = def_override.(fn{1});
    end
end

% Library argument setting the parameter nam to fraction f of its limit lim
par = @(nam,lim,f) {sprintf('%s=%.6g;',nam,f*lim)};

% Settings for V(x) = <Zop*Z(x),Pmat*Zop*Z(u)>, smallest degrees in piesos_ineq
% that work: (d,p)=(1,2), see STATUS
s_pie = struct('lf','pie','ineq_degs',struct('in',1,'out',2));

% Examples: {number, expected, settings, library arguments, fraction}
ex_list = 34;
ex_tab = {
    33,  true,  struct(), {},                       NaN;    % transport, c1=c2=-1
    33,  true,  struct(), {'c1=-2;c2=-0.5;'},       NaN;
    34,  true,  struct(), par('lam',2*pi^2,0.75),   0.75;
    34,  true,  struct(), par('lam',2*pi^2,0.85),   0.85;
    34,  NaN,   struct(), par('lam',2*pi^2,0.90),   0.90;
    34,  NaN,   struct(), par('lam',2*pi^2,0.95),   0.95;
    34,  false, struct(), par('lam',2*pi^2,1.05),   1.05;
    35,  true,  struct(), {},                       NaN;    % advection-diffusion
    35,  true,  struct(), {'b1=2;b2=4;'},           NaN;
    37,  true,  struct(), par('a',pi^2/2,0.50),     0.50;   % Dirichlet-Neumann
    37,  true,  struct(), par('a',pi^2/2,0.75),     0.75;
    37,  true,  struct(), par('a',pi^2/2,0.99),     0.99;
    37,  false, struct(), par('a',pi^2/2,1.05),     1.05;
    36,  true,  s_pie,    {},                       NaN;    % telegraph, damped
    36,  false, s_pie,    {'lam=-1;'},              NaN;    % negative damping
    38,  true,  s_pie,    {},                       NaN;    % wave, undamped
    40,  true,  struct(), {},                       NaN;    % ODE-PDE, K=-2
    40,  true,  struct(), {'K=-1.25;'},             NaN;
    40,  false, struct(), {'K=-0.5;'},              NaN;    % ODE unstable
    42,  true,  struct(), {},                       NaN;    % ODE, heat, delay
    43,  NaN,   struct(), {},                       NaN;    % +a1*x(t-tau)
    43,  true,  struct(), {'a1=-1;'},               NaN;    % -a1*x(t-tau)
    };

% Rows of the table to run
if ~exist('ex_list','var')
    ex_list = unique(cell2mat(ex_tab(:,1)))';
end
rows = find(any(abs(cell2mat(ex_tab(:,1))-ex_list(:)')<1e-9,2))';
if isempty(rows)
    error("None of the examples in ex_list are in the table ex_tab.")
end




%%%% 1. Run the examples

n_ex = numel(rows);
res = struct('label',cell(1,n_ex),'status',[],'feasratio',NaN,'time',[],'peak',NaN,'sdp','','msg',[]);
for k = 1:n_ex
    row = rows(k);
    ex_num = ex_tab{row,1};
    opt = def;
    for fn = fieldnames(ex_tab{row,3})'
        opt.(fn{1}) = ex_tab{row,3}.(fn{1});
    end
    if ~ismember(opt.lf,{'pde','pie'})
        error("Setting lf must be 'pde' or 'pie'.")
    end
    res(k).label = strtrim(sprintf('%g %s',ex_num,strjoin(ex_tab{row,4})));
    if strcmp(opt.lf,'pie')
        res(k).label = [res(k).label ' [pie]'];
    end
    fprintf('\n===== Example %s (%d/%d) =====\n',res(k).label,k,n_ex);
    t0 = tic;
    try
        % Declare the PIE
        PDE = examples_PDE_library_PIETOOLS_nouinput(ex_num,ex_tab{row,4}{:});
        PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');

        % Options for piesos_ineq
        ineq_opts = struct('psatz',opt.psatz);
        if ~isempty(opt.ineq_degs)
            ineq_opts.degs = opt.ineq_degs;
        end

        % Declare V and its derivative dV
        x = PIE.x;
        u = PIE.T*x;
        prog = piesos_program(x);
        Popts.degs = opt.pdegs;
        Popts.exclude = opt.exclude;
        if strcmp(opt.lf,'pde')
            % V(x) = <Zop*Z(u),Pmat*Zop*Z(u)>, with Pmat>=0
            Popts.pos = true;
            [prog,Vx,dV] = piesos_LF(prog,PIE,'pde',(1:opt.d)',Popts);

            % Add positivity term eppos*<u,u>, with derivative 2*eppos*<u,f(x)>
            Vx = Vx + opt.eppos*innerprod(u,u);
            dV = dV + 2*opt.eppos*innerprod(u,PIE.f);
        else
            % V(x) = <Zop*Z(x),Pmat*Zop*Z(u)>, with symmetric Pmat. Positivity
            % V>=eppos*<u,u> is enforced here, as in Cor. 35 of Jagt & Peet (2025)
            Popts.pos = false;
            [prog,Vx,dV] = piesos_LF(prog,PIE,'pie',(1:opt.d)',Popts);
            prog = piesos_ineq(prog,Vx-opt.eppos*innerprod(u,u),ineq_opts);
        end

        % Enforce dV <= 0 and solve
        prog = piesos_ineq(prog,-dV,ineq_opts);
        if isempty(opt.solver)
            prog = piesos_solve(prog);
        else
            prog = piesos_solve(prog,struct('solver',opt.solver));
        end

        info = prog.solinfo.info;
        res(k).feasratio = info.feasratio;
        res(k).sdp = sprintf('%dx%d',numel(prog.solinfo.y),numel(prog.solinfo.x));
        if info.pinf==0 && info.dinf==0 && info.feasratio>=opt.feas_tol
            res(k).status = 'feasible';
        else
            res(k).status = 'infeasible';
        end
        res(k).msg = sprintf('pinf=%d, dinf=%d, numerr=%d',info.pinf,info.dinf,info.numerr);
    catch ME
        res(k).status = 'error';
        res(k).msg = ME.message;
        fprintf(2,'Example %s failed: %s\n',res(k).label,ME.message);
        for i_st = 1:min(4,numel(ME.stack))
            fprintf(2,'    in %s (line %d)\n',ME.stack(i_st).name,ME.stack(i_st).line);
        end
    end
    res(k).time = toc(t0);
    try
        % Peak memory of the Matlab process so far (Windows only)
        res(k).peak = double(System.Diagnostics.Process.GetCurrentProcess().PeakWorkingSet64)/2^30;
    catch
    end
end




%%%% 2. Report

n_fail = 0;
fprintf('\n%-30s %-11s %-9s %-12s %-9s %-9s %-10s %s\n','Example','Status','feasratio','SDP (m x n)','Time [s]','Peak [GB]','Expected','Verdict');
for k = 1:n_ex
    expected = ex_tab{rows(k),2};
    if strcmp(res(k).status,'error')
        verdict = 'FAIL (error)';
    elseif isnan(expected)
        verdict = 'n/a';
    elseif expected==strcmp(res(k).status,'feasible')
        verdict = 'ok';
    else
        verdict = 'FAIL';
    end
    n_fail = n_fail + startsWith(verdict,'FAIL');
    exp_str = 'unknown';
    if ~isnan(expected)
        exp_str = 'unstable';
        if expected
            exp_str = 'stable';
        end
    end
    fprintf('%-30s %-11s %-9.3g %-12s %-9.1f %-9.2f %-10s %s  %s\n',res(k).label,res(k).status, ...
            res(k).feasratio,res(k).sdp,res(k).time,res(k).peak,exp_str,verdict,res(k).msg);
end
fprintf('\n%d of %d examples failed.\n',n_fail,n_ex);

% Largest feasible fraction of the stability limit, per example with a sweep
% (peak is the peak of the whole Matlab process so far, not of one example)
fracs = cell2mat(ex_tab(rows,5));
for ex_num = unique(cell2mat(ex_tab(rows(~isnan(fracs)),1)))'
    in_sweep = abs(cell2mat(ex_tab(rows,1))-ex_num)<1e-9 & ~isnan(fracs);
    is_feas = strcmp({res.status},'feasible');
    fprintf('Example %g: largest feasible fraction of the stability limit is %g\n', ...
            ex_num,max([0,fracs(in_sweep & is_feas(:))']));
end
