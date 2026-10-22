%
% Test the SOS stability test on the 1D linear PDEs in the PDE library
%
% This script loops over the 1D stability examples of the library
% `examples_PDE_library_PIETOOLS_nouinput`, and for each one repeats the
% steps of `Test_Linear_Examples`: it converts the PDE to a PIE, declares
% V(x) = <Zop*Z(u),Pmat*Zop*Z(u)> + eppos*||u||^2 with `piesos_LF`, and
% enforces dV <= 0 with `piesos_ineq`. Each program is solved, and its
% status is compared to the expected stability of the example.
%
% Examples included: all 1D stability examples (1-19)
% Example 15 is unstable, and must be infeasible;
% for 6 and 7 the parameter lies just below the stability limit; 
% for 11 and 12 the parameter is beyond the bound proven in the literature,
% so the outcome is not prescribed. 
% For 4.5 the default sig2=2.048 is beyond the limit of 1.049 stated in the 
% library (probably a typo for 1.048), so it is also run with sig2=1.04. 
% For 3 (Diagne) the classical PIETOOLS stability test 
% (PIETOOLS_PIE2PDEstability, 'heavy') is infeasible as well, so the 
% outcome is not prescribed; 
% the same holds for 8 (ODE-PDE, Das).
%
% The user may change the list ex_list of examples to run (set it before
% running the script), the default settings `def` (also by a struct
% def_override), or a row of table `ex_tab`, with columns
%   {number, expected, settings, arguments},
% where expected is true (stable, must be feasible), false (unstable,
% must be infeasible) or NaN (not prescribed), settings is a struct
% with any of the fields of `def`, overriding the defaults, and arguments
% a cell of extra inputs for the library, e.g. {'sig2=1.04;'}.
%
% STATUS (10/22/2026, MOSEK, default settings, 29 rows): all rows give the
% expected result (feasible for the stable examples, infeasible for 4.5 with
% sig2=2.048, 15 and 16, no verdict for 3, 8, 11, 12), except 18. Notes:
% - 10 (ODE-PDE, Tang) is feasible with the LF in the PIE state, lf='pie'
%   (feasratio 1.0, ineq degs 3; also 0.96-1.01 for ineq degs 2-3, psatz
%   0-1), as in the classical test (PIETOOLS_PIE2PDEstability, feasible for
%   'light', 'heavy', 'veryheavy'). With lf='pde' it is infeasible (numerr=2
%   for k=-1 to -16 and all variations tried); V and dV are correct for this
%   PIE (Test_piesos_LF, case_ode_pde);
% - 16 (4th order Timoshenko beam) is unstable for the default gam=1: the
%   symbol of w_tttt+alp*w_ssss+bet*w_tt-gam*w_ttss, mu^2+(gam*k^2+bet)*mu
%   +alp*k^4=0 for mu=lam^2, has complex mu for large k if gam^2<4*alp, so
%   Re(lam) grows like k^2. The library also codes BC 2 as w(1)=0 instead of
%   the documented w_s(0)=0 (PIETOOLS_PDE_Ex_Timoshenko_Beam_3). The result
%   is infeasible for the new recipe (pdegs 1-3) and for the classical
%   primal and dual tests ('light', 'heavy': pinf=1; 'veryheavy': numerr=2),
%   also for gam=2 (physical value for unit parameters) with either BC 2,
%   whose stability is not verified;
% - 18 (Datko wave equation, stable by its characteristic equation
%   (1-z)exp(z)=1, z=2*(lam+a)): infeasible for the primal tests, new recipe
%   with lf='pde' (pdegs 1-3, ineq degs 2-3, psatz 0-2) and lf='pie', and the
%   classical primal test; the classical DUAL test is feasible (feasratio
%   0.98-1.08 for 'light', 'heavy', 'veryheavy'). Dual LFs are not available
%   in the new tools. Example 19, with u_s as a state, is feasible.
% Not supported: Lyapunov functionals of degree larger than 1 in
% vector-valued states (error in intop/combine_terms).
%
% DJ (CC), 10/22/2026: Initial coding;

%%%% 0. Settings

% Defaults for the Lyapunov functional and the inequality, see
% Test_Linear_Examples
def.lf = 'pde';         % LF in the PDE state ('pde', V(x)=<Z(Tx)P Z(Tx)>), or in the PIE state x ('pie', V(x)=<Z(x),PZ(Tx)>)
def.d = 1;              % degree of the monomial basis Z(u), doubled in V
def.pdegs = 2;          % maximal monomial degree of Zop
def.exclude = [0;0;0];  % operators of Zop to exclude (mult, lower, upper)
def.eppos = 1e-4;       % coefficient of the positivity term eppos*||u||^2
def.psatz = 1;          % order of the multipliers in piesos_ineq
def.ineq_degs = [];     % degrees of the monomials in piesos_ineq (default if [])
def.solver = '';        % SDP solver passed to piesos_solve ('' for the default)
def.feas_tol = 0.5;     % minimal feasratio to accept a solution
if exist('def_override','var')
    % Optionally override the defaults by a struct def_override
    for fn = fieldnames(def_override)'
        def.(fn{1}) = def_override.(fn{1});
    end
end

% Examples: {number, expected, settings, library arguments}
ex_list = 6;    % list of examples to run
ex_tab = {
    1,   true,  struct(), {};               % transport
    2.0, true,  struct(), {};               % Lamare, 2x2 balance law
    2.1, true,  struct(), {};
    2.2, true,  struct(), {};
    2.3, true,  struct(), {};
    2.4, true,  struct(), {};
    3,   NaN,   struct(), {};               % Diagne, Saint-Venant-Exner
    4.0, true,  struct(), {};               % Saba, 2x2 hyperbolic
    4.1, true,  struct(), {};
    4.2, true,  struct(), {};
    4.3, true,  struct(), {};
    4.4, true,  struct(), {};
    4.5, false, struct(), {};               % sig2=2.048 > 1.049
    4.5, true,  struct(), {'sig2=1.04;'};
    5,   true,  struct(), {};               % reaction-diffusion, lam<pi^2
    6,   true,  struct(), {};               % lam<pi^2/4
    7,   true,  struct(), {};               % variable coefficients
    8,   NaN,   struct(), {};               % ODE-PDE, Das
    9,   true,  struct(), {};               % ODE in BC
    10,  true,  struct('lf','pie','ineq_degs',struct('in',3,'out',3)), {};   % ODE-PDE, Tang
    11,  NaN,   struct(), {};               % Ahmadi, R=2.9
    12,  NaN,   struct(), {};               % Ahmadi, R=30
    13,  true,  struct(), {};               % Euler-Bernoulli beam
    14,  true,  struct(), {};               % Timoshenko, hyperbolic
    15,  false, struct(), {};               % Timoshenko, unstable
    16,  false, struct(), {};               % Timoshenko, 4th order, gam=1
    17,  true,  struct(), {};               % wave, boundary damped
    18,  true,  struct(), {};               % Datko
    19,  true,  struct(), {};
    };

% Rows of the table to run
if ~exist('ex_list','var')
    ex_list = cell2mat(ex_tab(:,1))';
end
rows = find(any(abs(cell2mat(ex_tab(:,1))-ex_list(:)')<1e-9,2))';
if isempty(rows)
    error("None of the examples in ex_list are in the table ex_tab.")
end




%%%% 1. Run the examples

n_ex = numel(rows);
res = struct('label',cell(1,n_ex),'status',[],'feasratio',NaN,'time',[],'msg',[]);
for k = 1:n_ex
    row = rows(k);
    ex_num = ex_tab{row,1};
    opt = def;
    for fn = fieldnames(ex_tab{row,3})'
        opt.(fn{1}) = ex_tab{row,3}.(fn{1});
    end
    res(k).label = strtrim(sprintf('%g %s',ex_num,strjoin(ex_tab{row,4})));
    fprintf('\n===== Example %s (%d/%d) =====\n',res(k).label,k,n_ex);
    t0 = tic;
    try
        % Declare the PIE
        PDE = examples_PDE_library_PIETOOLS_nouinput(ex_num,ex_tab{row,4}{:});
        PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');

        % Declare V and its derivative dV
        x = PIE.x;
        u = PIE.T*x;
        prog = piesos_program(x);
        Popts = struct('degs',opt.pdegs,'exclude',opt.exclude);
        ineq_opts = struct('psatz',opt.psatz);
        if ~isempty(opt.ineq_degs)
            ineq_opts.degs = opt.ineq_degs;
        end
        if strcmp(opt.lf,'pde')
            % V(u) = <Zop*Z(u),Pmat*Zop*Z(u)>, Pmat>=0, plus eppos*<u,u>,
            % with derivative 2*eppos*<u,f(x)>
            Popts.pos = true;
            [prog,Vx,dV] = piesos_LF(prog,PIE,'pde',(1:opt.d)',Popts);
            dV = dV + 2*opt.eppos*innerprod(u,PIE.f);
        else
            % V(x) = <Zop*Z(x),Pmat*Zop*Z(u)>, symmetric Pmat, as in the
            % classical test (V = <x,P*x>, P = Q'*T >= eppos*T'*T): impose
            % V >= eppos*<u,u> explicitly
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
end




%%%% 2. Report

n_fail = 0;
fprintf('\n%-18s %-11s %-9s %-9s %-10s %s\n','Example','Status','feasratio','Time [s]','Expected','Verdict');
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
    fprintf('%-18s %-11s %-9.3g %-9.1f %-10s %s  %s\n',res(k).label,res(k).status, ...
            res(k).feasratio,res(k).time,exp_str,verdict,res(k).msg);
end
fprintf('\n%d of %d examples failed.\n',n_fail,n_ex);
