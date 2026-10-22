%
% Test the SOS stability test on the DDE stability examples in the library
%
% This script runs the stability examples of `examples_DDE_library_PIETOOLS`
%   x'(t) = A0*x(t) + sum_j Ai{j}*x(t-tau(j))
%                   + sum_j Adi{j}*int_{-tau(j)}^{0} x(t+theta) dtheta,
% which are only given as commented code in that file, and are therefore
% repeated in the table `ex_tab` below, with the parameters of the library.
% Each DDE is written as an ODE-PDE with transport equations for the
% delayed states z_j(t,s) = x(t-tau(j)*s), s in [0,1],
%   x'      = A0*x + sum_j Ai{j}*z_j(1) + tau(j)*Adi{j}*int_0^1 z_j(s) ds;
%   z_j'    = -(1/tau(j))*z_j_{s};      z_j(0) = x;
% see also example 44 of `examples_PDE_library_PIETOOLS_nouinput`. The
% ODE-PDE is converted to a PIE, and stability is tested with the same SOS
% Lyapunov functional as in `Test_Linear_Examples_1D`.
%
% Each example is run at the delay of the library, which is stable, and
% must be feasible. Examples marked 'beyond' are run at a delay just beyond
% the stability limit reported in the library; these are unstable, so the
% program must be infeasible (a certificate would be a bug).
%
% The user may change the list ex_list of rows of ex_tab to run (set it
% before running the script), the defaults `def`, or the settings in a row
% of ex_tab, with columns {name, expected, DDE, settings}, where expected
% is true (stable), false (unstable) or NaN (not prescribed).
%
% STATUS (10/21/2026, MOSEK, default settings): all 14 rows give the expected
% result: the stable delays are feasible (feasratio >= 0.99) and the delays
% beyond the stability limit are infeasible. This needed support for vector-
% valued and ODE states in the SOS tools (see the log lines of 10/21/2026 in
% piesos_LF, piesos_ineq, @tensopvar/subsref and mtimes).
%
% DJ (CC), 10/22/2026: Initial coding;

%%%% 0. Settings

def.d = 1;              % degree of the monomial basis Z(u), doubled in V
def.pdegs = 2;          % maximal monomial degree of Zop
def.exclude = [0;0;0];  % operators of Zop to exclude (mult, lower, upper)
def.eppos = 1e-4;       % coefficient of the positivity term eppos*||u||^2
def.psatz = 1;          % order of the multipliers in piesos_ineq
def.ineq_degs = [];     % degrees of the monomials in piesos_ineq (default if [])
def.feas_tol = 0.5;     % minimal feasratio to accept a solution
if exist('def_override','var')
    % Optionally override the defaults by a struct def_override
    for fn = fieldnames(def_override)'
        def.(fn{1}) = def_override.(fn{1});
    end
end

ex_tab = {};
% 1. Distributed delay; stable for tau in [0.2, 2.04]
D = struct('A0',[0.2 0;0.2 0.1],'Ai',{{[0 0;0 0]}},'Adi',{{[-1 0;-1 -1]}},'tau',2.0);
ex_tab(end+1,:) = {'dist. delay, tau=2.0',true,D,struct()};
D.tau = 2.1;
ex_tab(end+1,:) = {'dist. delay, tau=2.1 (beyond)',false,D,struct()};
% 2. Distributed delay [1, p.270]; stable for tauN < 2
Ap1 = [2 2.25;0 -.5];   Ap2 = [-1 0;0 -1];   tauN = 1.8;
D = struct('A0',[-1.5 0;0.5 -1],'Adi',{{Ap1-Ap2,Ap2}},'tau',[tauN/2,tauN]);
ex_tab(end+1,:) = {'dist. delay [1], tauN=1.8',true,D,struct()};
% 3. Example B [1, p.193]; stable for tau in [0.100174, 1.71785]
D = struct('A0',[0 1;-2 .1],'Ai',{{[0 0;1 0]}},'tau',1.5);
ex_tab(end+1,:) = {'Example B, tau=1.5',true,D,struct()};
D.tau = 1.8;
ex_tab(end+1,:) = {'Example B, tau=1.8 (beyond)',false,D,struct()};
% 4. Example [1, p.192]; stable for tau < 6.17258
D = struct('A0',[-2 0;0 -.9],'Ai',{{[-1 0;-1 -1]}},'tau',6);
ex_tab(end+1,:) = {'[1, p.192], tau=6',true,D,struct()};
D.tau = 6.5;
ex_tab(end+1,:) = {'[1, p.192], tau=6.5 (beyond)',false,D,struct()};
% 5. Example A [2]; x' = -x(t-tau), stable for tau < 1.5707
D = struct('A0',0,'Ai',{{-1}},'tau',1.56);
ex_tab(end+1,:) = {'Example A, tau=1.56',true,D,struct()};
D.tau = 1.65;
ex_tab(end+1,:) = {'Example A, tau=1.65 (beyond)',false,D,struct()};
% 6. Example C [2,3]; x' = -2x + b*x(t-1) - x(t-2), stable for b < 3
D = struct('A0',-2,'Ai',{{2.99,-1}},'tau',[1,2]);
ex_tab(end+1,:) = {'Example C, b=2.99',true,D,struct()};
D.Ai{1} = 3.1;
ex_tab(end+1,:) = {'Example C, b=3.1 (beyond)',false,D,struct()};
% 7. Example D [2]; stable for tauN < 1.3722
tauN = 1.3;
D = struct('A0',[0 1;-1 .1],'Ai',{{[0 0;-1 0],[0 0;1 0]}},'tau',[tauN/2,tauN]);
ex_tab(end+1,:) = {'Example D, tauN=1.3',true,D,struct()};
tauN = 1.45;
D.tau = [tauN/2,tauN];
ex_tab(end+1,:) = {'Example D, tauN=1.45 (beyond)',false,D,struct()};
% 8. Example E [2,9]; stable for tau up to ~3.9 (needs 'veryheavy' in
% the old toolbox), run at the library delay 1.1
A = [0 0 1 0; 0 0 0 1;-10 10 0 0;5 -15 0 -.25];   B = [0;0;1;0];
C = [1 0 0 0];   K = 1;
D = struct('A0',A-B*K*C,'Ai',{{B*K*C}},'tau',1.1);
ex_tab(end+1,:) = {'Example E, tau=1.1',true,D,struct()};

if ~exist('ex_list','var')
    ex_list = 1:size(ex_tab,1);
end




%%%% 1. Run the examples

pvar s
n_ex = numel(ex_list);
res = struct('status',cell(1,n_ex),'feasratio',NaN,'time',[],'msg',[]);
for k = 1:n_ex
    row = ex_list(k);
    DDE = ex_tab{row,3};
    opt = def;
    for fn = fieldnames(ex_tab{row,4})'
        opt.(fn{1}) = ex_tab{row,4}.(fn{1});
    end
    fprintf('\n===== DDE %d: %s (%d/%d) =====\n',row,ex_tab{row,1},k,n_ex);
    t0 = tic;
    try
        % Declare the ODE-PDE with transport equations for the delays
        n = size(DDE.A0,1);
        nd = numel(DDE.tau);
        xo = pde_var(n);
        z = cell(nd,1);
        for j = 1:nd
            z{j} = pde_var(n,s,[0,1]);
        end
        rhs = DDE.A0*xo;
        for j = 1:nd
            if isfield(DDE,'Ai') && numel(DDE.Ai)>=j && ~isempty(DDE.Ai{j})
                rhs = rhs + DDE.Ai{j}*subs(z{j},s,1);
            end
            if isfield(DDE,'Adi') && numel(DDE.Adi)>=j && ~isempty(DDE.Adi{j})
                rhs = rhs + DDE.tau(j)*DDE.Adi{j}*int(z{j},s,[0,1]);
            end
        end
        PDE = diff(xo,'t')==rhs;
        for j = 1:nd
            PDE = [PDE; diff(z{j},'t')==-(1/DDE.tau(j))*diff(z{j},s);
                        subs(z{j},s,0)==xo];
        end
        PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');

        % Declare V and its derivative dV
        x = PIE.x;
        prog = piesos_program(x);
        Popts.degs = opt.pdegs;
        Popts.exclude = opt.exclude;
        Popts.pos = true;
        [prog,Vx,dV] = piesos_LF(prog,PIE,'pde',(1:opt.d)',Popts);

        % Add positivity term eppos*<u,u>, with derivative 2*eppos*<u,f(x)>
        u = PIE.T*x;
        Vx = Vx + opt.eppos*innerprod(u,u);
        dV = dV + 2*opt.eppos*innerprod(u,PIE.f);

        % Enforce dV <= 0 and solve
        ineq_opts = struct('psatz',opt.psatz);
        if ~isempty(opt.ineq_degs)
            ineq_opts.degs = opt.ineq_degs;
        end
        prog = piesos_ineq(prog,-dV,ineq_opts);
        prog = piesos_solve(prog);

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
        fprintf(2,'DDE %d failed: %s\n',row,ME.message);
        for i_st = 1:min(4,numel(ME.stack))
            fprintf(2,'    in %s (line %d)\n',ME.stack(i_st).name,ME.stack(i_st).line);
        end
    end
    res(k).time = toc(t0);
end




%%%% 2. Report

n_fail = 0;
fprintf('\n%-34s %-11s %-9s %-9s %-10s %s\n','Example','Status','feasratio','Time [s]','Expected','Verdict');
for k = 1:n_ex
    row = ex_list(k);
    expected = ex_tab{row,2};
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
    fprintf('%-34s %-11s %-9.3g %-9.1f %-10s %s  %s\n',ex_tab{row,1},res(k).status, ...
            res(k).feasratio,res(k).time,exp_str,verdict,res(k).msg);
end
fprintf('\n%d of %d examples failed.\n',n_fail,n_ex);
