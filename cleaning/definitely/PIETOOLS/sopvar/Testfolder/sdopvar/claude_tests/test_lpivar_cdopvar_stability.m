%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_LPIVAR_CDOPVAR_STABILITY runs the Q-FORM stability LPI of
% 'PIETOOLS_PIE2PDEstability' - the default 'stability' executive - through
% the container path, and checks it against the stock opvar path on the
% same PIE. It is the end-to-end test for 'lpivar_cdopvar':
%
%   poscopvar           the Lyapunov operator P
%   lpivar_cdopvar      the general (non-self-adjoint) operator Q
%   @cdopvar/minus      Top'*Q - P, constrained to zero WITHOUT 'symmetric'
%   @cdopvar/mtimes     A'*Q + Q'*A, a sum of fixed x decision products
%   @cdopvar/plus       with the negativity slack from poscopvar
%   lpi_eq_cdopvar      both equality constraints
%
% against the stock sequence with 'poslpivar', 'lpivar' and 'lpi_eq'.
% Degrees: P and the slack are given in poslpivar's vocabulary and mapped
% across as in 'test_poscopvar_stability'; Q's degrees come from the stock
% 'get_lpivar_degs(Pop,Top)' and are handed unchanged to both paths, which
% is lpivar_cdopvar's legacy [d1 d2 d3] mapping. The two paths therefore
% pose the same semidefinite problem; their decision variable counts must be
% equal and their certified verdicts must agree.
%
% Verdicts are GRADED. +1 is feasible at numerr = 0, and also at numerr = 1
% with |feasratio - 1| < 0.1, printed 'feasible*'; -1 is infeasible (pinf = 1)
% at numerr = 0, and also at numerr = 1 with |feasratio + 1| < 0.1, printed
% 'infeasible*'; anything else is 0, uncertified. SeDuMi's numerr = 1 is
% "desired accuracy not reached", not failure. The stock Q-form returns it on
% every stable system at every eppos2 in [1e-6, 1e-2] and both degree sets
% tried (feasratio +0.99 to +1.02), so requiring numerr = 0 would make the
% test vacuous; and the two paths, posing the same problem with differently
% ordered data, can land on either side of the accuracy threshold for the
% same verdict (measured: unstable PDE, stock numerr 0 vs container numerr 1,
% both pinf = 1, feasratio -1.02 / -1.03). This is looser than
% 'test_copvar_kyp', where numerr = 1 came with feasratio anywhere in
% [-0.4, +0.7] and certified nothing. Stable systems must be +1 and unstable
% ones -1 on both paths; the one case with no known answer only has to agree.
%
% Systems, on x(t) in R and u(s,t) in L_2[0,1], from test_poscopvar_stability:
%   xdot = -c x + kappa int_0^1 u ds,   u_t = u_ss + lam u - kappa x,
%   u(0) = u(1) = 0; energy-neutral coupling, stable for c > 0, lam = 0.
%
% MMP, 09/25/2026: Initial coding
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');

sopts = struct('solver','sedumi');
eppos2 = 1e-6;
dd1 = {1,[1,1,2],[1,1,2]};
dd2 = {2,[2,2,4],[2,2,4]};

cases = {
    {'c=1  lam=0   kappa=1  stable, coupled',      1,  0,  1, +1}
    {'c=1  lam=0   kappa=0  stable, decoupled',    1,  0,  0, +1}
    {'c=-1 lam=0   kappa=0  unstable ODE',        -1,  0,  0, -1}
    {'c=1  lam=12  kappa=0  unstable PDE',         1, 12,  0, -1}
    {'c=1  lam=12  kappa=1  unstable, coupled',    1, 12,  1,  0}
    };

npass = 0;
for ic = 1:numel(cases)
    [lbl,c,lam,kap,want] = deal(cases{ic}{:});
    PIE = build_pie(c,lam,kap);
    [s_op,i_op,n_op,Qdeg] = run_opvar(PIE,dd1,dd2,eppos2,sopts);
    [s_co,i_co,n_co] = run_copvar(PIE,dd1,dd2,eppos2,sopts,Qdeg);
    if n_op~=n_co
        error(['test_lpivar_cdopvar_stability: ''%s'': the programs declare %d and %d '...
               'decision variables; the families should coincide.'],lbl,n_op,n_co);
    end
    if s_op~=s_co
        error(['test_lpivar_cdopvar_stability: ''%s'': stock path %s, container path %s.'], ...
              lbl,i_op,i_co);
    end
    if want~=0 && s_op~=want
        error(['test_lpivar_cdopvar_stability: ''%s'' should be certified %s, but both '...
               'paths report %s.'],lbl,ternary(want>0,'feasible','infeasible'),i_op);
    end
    fprintf('  passed: %-38s %-11s Qdeg %s  %d dvars  (opvar %s | copvar %s)\n',lbl, ...
            verdict(s_op,i_op),mat2str(full(Qdeg)),n_op,i_op,i_co);
    npass = npass+1;
end
fprintf('lpivar_cdopvar stability test passed (%d of %d cases).\n',npass,numel(cases));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE = build_pie(c,lam,kap)
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
function [st,info,nd,Qdeg] = run_opvar(PIE,dd1,dd2,eppos2,sopts)
% 'PIETOOLS_PIE2PDEstability' with override = 1, psatz = 0, epneg = 0 and
% the equality form of the negativity constraint.
Top = PIE.T;    Aop = PIE.A;
opt0 = struct('psatz',0,'exclude',[0,0,0,0],'sep',0);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pop] = poslpivar(prog,Top.dim,dd1,opt0);
Pop = Pop + eppos2*Top'*Top;
Qdeg = get_lpivar_degs(Pop,Top);
[prog,Qop] = lpivar(prog,Top.dim,Qdeg);
prog = lpi_eq(prog,Top'*Qop-Pop);
Dop = Aop'*Qop + Qop'*Aop;
[prog,Deop] = poslpivar(prog,Dop.dim,dd2,opt0);
prog = lpi_eq(prog,Dop+Deop,'symmetric');
nd = numel(prog.decvartable);
evalc('sol = lpisolve(prog,sopts);');
[st,info] = read_sol(sol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st,info,nd] = run_copvar(PIE,dd1,dd2,eppos2,sopts,Qdeg)
% The same LPI through the containers.
Tm = opvar2copvar(PIE.T);
Am = opvar2copvar(PIE.A);
n = PIE.T.dim(1,1);     m = PIE.T.dim(2,1);
vname = PIE.T.var1.varname;
spaces = {{},reshape(vname,1,[])};      % R^n and L_2^m[s]
dims = [n;m];

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pm] = poscopvar(prog,dims,spaces,PIE.dom,pl2pm(dd1));
Pm = Pm + eppos2*(Tm'*Tm);
[prog,Qm] = lpivar_cdopvar(prog,dims,spaces,PIE.dom,Qdeg);
prog = lpi_eq_cdopvar(prog,Tm'*Qm - Pm);               % NOT symmetric
Dm = (Am')*Qm + (Qm')*Am;
[prog,Nm] = poscopvar(prog,dims,spaces,PIE.dom,pl2pm(dd2));
prog = lpi_eq_cdopvar(prog,Dm+Nm,'symmetric');
nd = numel(prog.decvartable);
evalc('sol = lpisolve(prog,sopts);');
[st,info] = read_sol(sol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = pl2pm(dd)
% 'poslpivar' degrees to 'poscopvar' degrees, as in test_poscopvar_stability.
d1 = dd{1};     d2 = dd{2};     d3 = dd{3};
deg = { struct('int',0), ...
        { struct('int',d1,'mult',0), ...
          struct('int',d2(1),'mult',d2(2),'joint',d2(3)), ...
          struct('int',d3(1),'mult',d3(2),'joint',d3(3)) } };
end


function [st,info] = read_sol(sol)
% Graded verdict, see the header: +1 feasible (numerr 0, or numerr 1 with
% feasratio within 0.1 of +1), -1 certified infeasible, 0 otherwise.
pinf = NaN;     ne = NaN;       fr = NaN;
if isfield(sol,'solinfo') && isfield(sol.solinfo,'info')
    q = sol.solinfo.info;
    if isfield(q,'pinf'),       pinf = q.pinf;      end
    if isfield(q,'numerr'),     ne = q.numerr;      end
    if isfield(q,'feasratio'),  fr = q.feasratio;   end
end
st = 0;
if pinf==0 && (ne==0 || (ne==1 && abs(fr-1)<0.1)),  st = +1;    end
if pinf==1 && (ne==0 || (ne==1 && abs(fr+1)<0.1)),  st = -1;    end
info = sprintf('pinf=%g numerr=%g fr=%+.2f',pinf,ne,fr);
end


function s = verdict(st,info)
switch st
    case +1
        s = 'feasible';
        if nargin>1 && contains(info,'numerr=1'),   s = 'feasible*';    end
    case -1
        s = 'infeasible';
        if nargin>1 && contains(info,'numerr=1'),   s = 'infeasible*';  end
    otherwise,  s = 'uncertified';
end
end


function s = ternary(c,a,b)
if c,   s = a;  else,   s = b;  end
end
