%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_POSCOPVAR_STABILITY runs a stability LPI on a MIXED domain through the
% container path, and checks it against the stock opvar path on the same PIE.
%
% This is the end-to-end test for 'poscopvar': a PIE whose state lives in
% R^n x L_2^m[a,b] cannot be handled by 'possopvar' at all, since that builds
% only the L_2 -> L_2 block. The workflow exercised here is
%
%   opvar2copvar        T and A into the container
%   poscopvar           the Lyapunov operator and the negativity operator
%   @cdopvar/mtimes     T'*P*A + A'*P*T on the mixed domain
%   @cdopvar/plus       the strict-positivity identity, and D + N
%   lpi_eq_cdopvar      the equality constraint
%
% against 'PIETOOLS_PDEstability's own sequence with 'poslpivar', 'lpi_eq'
% and the 4-PI 'opvar'. The degrees are given once in poslpivar's vocabulary
% and mapped across by the correspondence that
% 'test_poscopvar_vs_poslpivar' verifies, so the two runs pose the SAME
% semidefinite feasibility problem and must return the same verdict.
%
% System, on x(t) in R and u(s,t) in L_2[0,1]:
%
%   xdot = -c x + kappa int_0^1 u ds
%   u_t  = u_ss + lam u - kappa x,          u(0) = u(1) = 0
%
% The coupling is energy-neutral: with V = x^2 + int u^2,
%
%   Vdot = -2c x^2 - 2 int u_s^2 + 2 lam int u^2,
%
% so c>0 with lam=0 is stable for every kappa, while a decoupled case with
% c<0, or with lam>pi^2, is unstable. Those four are absolute anchors, so
% agreement between the two paths cannot be achieved by both being wrong in
% the same way. The fifth case has no independent verdict and is checked for
% agreement only.
%
% Requires SeDuMi. MOSEK reports these tests feasible when they are not - see
% the note in 'test_possopvar_integration' - so the solver is named
% explicitly rather than left to the default.
%
% MMP, 09/21/2026: Initial coding
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  lpi_eq_mdopvar -> lpi_eq_cdopvar,
%                  opvar2mopvar -> opvar2copvar, posmopvar -> poscopvar,
%                  run_mopvar -> run_copvar,
%                  test_posmopvar_stability -> test_poscopvar_stability,
%                  test_posmopvar_vs_poslpivar -> test_poscopvar_vs_poslpivar.
%                  File was 'test_posmopvar_stability.m'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

sopts = struct('solver','sedumi');
eppos = 1e-4;       eppos2 = 1e-6;

% Degrees, in 'poslpivar' vocabulary: d{1} for Z1 of the multiplier block,
% d{2} and d{3} for the two integral blocks.
dd1 = {1,[1,1,2],[1,1,2]};
dd2 = {2,[2,2,4],[2,2,4]};

% {label, c, lam, kappa, expected verdict or '' for agreement only}
cases = {
    {'c=1  lam=0   kappa=1  stable, coupled',      1,  0,  1, 'feas'  }
    {'c=1  lam=0   kappa=0  stable, decoupled',    1,  0,  0, 'feas'  }
    {'c=-1 lam=0   kappa=0  unstable ODE',        -1,  0,  0, 'infeas'}
    {'c=1  lam=12  kappa=0  unstable PDE',         1, 12,  0, 'infeas'}
    {'c=1  lam=12  kappa=1  unstable, coupled',    1, 12,  1, ''      }
    };

npass = 0;
for ic = 1:numel(cases)
    [lbl,c,lam,kap,want] = deal(cases{ic}{:});
    PIE = build_pie(c,lam,kap);

    [f_op,i_op] = run_opvar(PIE,dd1,dd2,eppos,eppos2,sopts);
    [f_mo,i_mo] = run_copvar(PIE,dd1,dd2,eppos,eppos2,sopts);

    if f_op~=f_mo
        error(['test_poscopvar_stability: the two paths disagree for case '...
               '''%s'': opvar says %s, copvar says %s.'],lbl,i_op,i_mo);
    end
    if ~isempty(want) && f_op~=strcmp(want,'feas')
        error(['test_poscopvar_stability: case ''%s'' should be %s, but both '...
               'paths report %s.'],lbl,want,i_op);
    end

    fprintf('  passed: %-38s %-8s (opvar %s, copvar %s)\n',lbl, ...
            ternary(f_op,'feasible','infeasible'),i_op,i_mo);
    npass = npass+1;
end

fprintf('poscopvar stability test passed (%d of %d cases).\n',npass,numel(cases));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE = build_pie(c,lam,kap)
% The coupled ODE-PDE system above, as a PIE.

pvar s1
X = pde_var();
u = pde_var(1,s1,[0,1]);
eqs = [diff(X,'t') == -c*X + kap*int(u,s1,[0,1]);
       diff(u,'t') == diff(u,s1,2) + lam*u - kap*X;
       subs(u,s1,0) == 0;
       subs(u,s1,1) == 0];
evalc('PIE = convert(eqs,''pie'');');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feas,info] = run_opvar(PIE,dd1,dd2,eppos,eppos2,sopts)
% The core of 'PIETOOLS_PDEstability': the equality form of the negativity
% constraint, one positive operator each (override = 1) and psatz = 0.

Top = PIE.T;    Aop = PIE.A;
opt0 = struct('psatz',0,'exclude',[0,0,0,0],'sep',0);

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pop] = poslpivar(prog,Top.dim,dd1,opt0);
Imat = blkdiag(eppos*eye(Pop.dim(1,:)),eppos2*eye(Pop.dim(2,:)));
Pop = Pop + mat2opvar(Imat,Pop.dim(:,2),PIE.vars,PIE.dom);
Dop = Top'*Pop*Aop + Aop'*Pop*Top;
[prog,Deop] = poslpivar(prog,Dop.dim,dd2,opt0);
prog = lpi_eq(prog,Dop+Deop,'symmetric');
evalc('sol = lpisolve(prog,sopts);');
[feas,info] = read_sol(sol);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feas,info] = run_copvar(PIE,dd1,dd2,eppos,eppos2,sopts)
% The same test through the container.

Tm = opvar2copvar(PIE.T);
Am = opvar2copvar(PIE.A);
n = PIE.T.dim(1,1);     m = PIE.T.dim(2,1);
vname = PIE.T.var1.varname;
spaces = {{},reshape(vname,1,[])};      % R^n and L_2^m[s]
dims = [n;m];

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = poscopvar(prog,dims,spaces,PIE.dom,pl2pm(dd1));
Imat = blkdiag(eppos*eye(n),eppos2*eye(m));
Pm = Pm + opvar2copvar(mat2opvar(Imat,[n;m],PIE.vars,PIE.dom));
Dm = (Tm')*(Pm*Am) + (Am')*(Pm*Tm);
[prog,Nm] = poscopvar(prog,dims,spaces,PIE.dom,pl2pm(dd2));
prog = lpi_eq_cdopvar(prog,Dm+Nm,'symmetric');
evalc('sol = lpisolve(prog,sopts);');
[feas,info] = read_sol(sol);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = pl2pm(dd)
% 'poslpivar' degrees to 'poscopvar' degrees. The R^n basis operator of
% poslpivar is the first row of Zop, which is the identity, hence degree 0 in
% the integration variable; over the L_2 space alpha=1 is Z1, alpha=2 is Z2
% and alpha=3 is Z3. This is the correspondence checked by
% 'test_poscopvar_vs_poslpivar', where the two families are shown to span the
% same set of operators.

d1 = dd{1};     d2 = dd{2};     d3 = dd{3};
deg = { struct('int',0), ...
        { struct('int',d1,'mult',0), ...
          struct('int',d2(1),'mult',d2(2),'joint',d2(3)), ...
          struct('int',d3(1),'mult',d3(2),'joint',d3(3)) } };

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [feas,info] = read_sol(sol)
% Primal feasibility as the solver reports it, with the numerical error flag
% alongside: a verdict reached with numerr>0 is not one to compare.

feas = true;    pinf = NaN;     ne = NaN;
if isfield(sol,'solinfo') && isfield(sol.solinfo,'info')
    q = sol.solinfo.info;
    if isfield(q,'pinf'),   pinf = q.pinf;      feas = ~(pinf==1);  end
    if isfield(q,'numerr'), ne = q.numerr;                          end
end
info = sprintf('pinf%d nerr%d',pinf,ne);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function o = ternary(c,a,b)

if c, o = a; else, o = b; end

end
