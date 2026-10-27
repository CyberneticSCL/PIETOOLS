%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST_COPVAR_KYP runs the H-infinity gain KYP LPI at FIXED gamma through
% the container path, and checks it against the stock opvar path on the
% same PIE. It is the end-to-end test for the container block operations:
%
%   opvar2copvar            T, A, B1, C1, D11 and the identities
%   poscopvar               the storage operator P and the negativity slack
%   @cdopvar/mtimes, plus   the KYP entries
%   @copvar/uminus          -(gam*Iw), -(gam*Iz)                            % MMP, 09/25/2026
%   @cdopvar/uminus         (-1)*N in the constraint: applied once, so      % MMP, 09/25/2026
%                           a sign error cannot cancel                      % MMP, 09/25/2026
% (was)   @cdopvar/uminus         -(gam*Iw), -(gam*Iz), and -KYP below      % MMP, 09/25/2026 (was)
%   @cdopvar/horzcat,vertcat  the 3 x 3 KYP operator, mixing fixed and
%                           decision blocks, R^n and L_2 spaces, and
%                           operators over different variable registries
%                           (Iw and Dzw have none)
%   @cdopvar/minus          KYP - ((-1)*N) = 0, i.e. KYP = -N <= 0          % MMP, 09/25/2026
% (was)   @cdopvar/minus          N - (-KYP) = 0, i.e. KYP = -N <= 0        % MMP, 09/25/2026 (was)
%   lpi_eq_cdopvar          the equality constraint
%
% against the construction of 'PIETOOLS_Hinf_gain_coercive' with 'opvar'
% arithmetic, 'poslpivar' and 'lpi_eq'. Both paths pose the SAME semidefinite
% feasibility problem - the storage and slack cones coincide, as
% 'test_poscopvar_vs_poslpivar' shows - so they must agree on every gamma.
%
% Only CERTIFIED verdicts are compared. This LPI is numerically hard near
% its threshold: measured on the ODE-PDE system, both paths return
% numerr=1 for every gamma in roughly [0.28, 1.2], with pinf verdicts that
% flip from point to point on either path, while below and above that band
% both are clean (numerr=0) and agree. A verdict with numerr=1 certifies
% nothing (the solver flags in 'test_poscopvar_stability' carry the same
% caveat), so a bisection on pinf alone lands inside the band. Instead the
% stock path is walked down by halving until it is cleanly infeasible twice
% in a row, and up by doubling until it is cleanly feasible twice; at those
% four gammas the container path must also be clean and agree. The band in
% between is reported, not asserted.
%
% Strictness is added identically on both paths: P >= eps I, and eps I on
% the finite-dimensional w and z blocks of the KYP operator, the latter
% built with the container 'blkdiag' (with a zero block on the state). NOT
% on the state block: strict negativity there, in the PIE's fundamental
% state norm, is infeasible for diffusion at every gamma (measured), which
% is why the executives leave it out.
%
% Structure is checked as well: the two programs must declare the same
% number of decision variables, since their cones coincide.
%
% Systems, on x(t) in R, u(s,t) in L_2[0,1], disturbance w, output z:
%   ODE-PDE     xdot = -x + int u,   u_t = u_ss - x + s w,   z = int u + x
%   PDE only    u_t = u_ss + s w,                            z = int u
%   feedthrough the ODE-PDE with z = int u + x + 0.5 w  (D11 ~= 0)
% all with u(0) = u(1) = 0. The KYP grid is 4 x 4 for the ODE-PDE cases
% (spaces w, z, R^1, L_2[s]) and 3 x 3 for the PDE-only case.
%
% MMP, 09/25/2026: Initial coding
% MMP, 09/25/2026: Constraint written KYP - ((-1)*N) instead of N - (-KYP).
%                  The old form negates KYP twice, so a @cdopvar/uminus that
%                  returned its input would still pass; now uminus acts once,
%                  on N. -(gam*Iw) is @copvar/uminus, not @cdopvar.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;
warning('off','sopvar:noncanonicalMultiplier');
warning('off','sdopvar:noncanonicalMultiplier');

sopts = struct('solver','sedumi');
% Degrees in 'poslpivar' vocabulary, mapped across by 'pl2pm'.
dd1 = {1,[1,1,2],[1,1,2]};
dd2 = {2,[2,2,4],[2,2,4]};

SYS = { 'ODE-PDE',      1, 0
        'PDE only',     0, 0
        'feedthrough',  1, 0.5 };

npass = 0;  nrun = 0;
for is = 1:size(SYS,1)
    [lbl,withode,d11] = deal(SYS{is,:});
    PIE = build_pie(withode,d11);

    % Clean anchors on the stock path: two certified-infeasible gammas going
    % down from 1 by halving, two certified-feasible ones going up by
    % doubling. 'st' is +1 clean feasible, -1 clean infeasible, 0 unclean.
    lows = [];      g = 1;
    while numel(lows)<2 && g>1e-4
        st = run_opvar(PIE,g,dd1,dd2,sopts);
        if st==-1,  lows(end+1) = g;    else,   lows = [];  end         %#ok<AGROW>
        g = g/2;
    end
    highs = [];     g = 1;
    while numel(highs)<2 && g<1e4
        st = run_opvar(PIE,g,dd1,dd2,sopts);
        if st==+1,  highs(end+1) = g;   else,   highs = []; end         %#ok<AGROW>
        g = 2*g;
    end
    if numel(lows)<2 || numel(highs)<2
        error('test_copvar_kyp: no clean anchors found for ''%s''.',lbl)
    end
    fprintf('  %-12s certified infeasible at %s, feasible at %s (stock path)\n',lbl, ...
            mat2str(lows,4),mat2str(highs,4));

    for g = [lows, highs]
        [s_op,i_op,n_op] = run_opvar(PIE,g,dd1,dd2,sopts);
        [s_co,i_co,n_co,K] = run_copvar(PIE,g,dd1,dd2,sopts);
        nrun = nrun+1;
        if s_co~=s_op
            error(['test_copvar_kyp: at gamma = %g for ''%s'' the stock path is %s but '...
                   'the container path is %s.'],g,lbl,i_op,i_co);
        end
        if n_co~=n_op
            error(['test_copvar_kyp: the two programs declare %d and %d decision '...
                   'variables for ''%s''; the cones should coincide.'],n_op,n_co,lbl);
        end
        fprintf('    gamma %-8.4g %-10s  KYP grid %s  %d dvars  (opvar %s | copvar %s)\n',g, ...
                ternary(s_op>0,'feasible','infeasible'),mat2str(size(K)),n_op,i_op,i_co);
        npass = npass+1;
    end
end

fprintf('test_copvar_kyp passed (%d of %d gamma points, %d systems).\n',npass,nrun,size(SYS,1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE = build_pie(withode,d11)
pvar s1
u = pde_var(1,s1,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
if withode
    X = pde_var();
    eqs = [diff(X,'t') == -X + int(u,s1,[0,1]);
           diff(u,'t') == diff(u,s1,2) - X + s1*w;
           z == int(u,s1,[0,1]) + X + d11*w;
           subs(u,s1,0) == 0;
           subs(u,s1,1) == 0];
else
    eqs = [diff(u,'t') == diff(u,s1,2) + s1*w;
           z == int(u,s1,[0,1]) + d11*w;
           subs(u,s1,0) == 0;
           subs(u,s1,1) == 0];
end
evalc('PIE = convert(eqs,''pie'');');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st,info,nd] = run_opvar(PIE,gam,dd1,dd2,sopts)
% 'PIETOOLS_Hinf_gain_coercive' at fixed gamma, Tw = 0, equality form, with
% the strictness described in the header.
ep = 1e-4;
Top = PIE.T;   Aop = PIE.A;   Bwop = PIE.B1;   Czop = PIE.C1;   Dzwop = PIE.D11;
opt0 = struct('psatz',0,'exclude',[0,0,0,0],'sep',0);
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[prog,Pop] = poslpivar(prog,Top.dim,dd1,opt0);
Pop = Pop + mat2opvar(ep*eye(sum(Top.dim(:,1))),Top.dim(:,2),PIE.vars,PIE.dom);
Iw = mat2opvar(eye(size(Bwop,2)),Bwop.dim(:,2),PIE.vars,PIE.dom);
Iz = mat2opvar(eye(size(Czop,1)),Czop.dim(:,1),PIE.vars,PIE.dom);
Dop = [-gam*Iw,           Dzwop',    (Pop*Bwop)'*Top;
        Dzwop,            -gam*Iz,   Czop;
        Top'*(Pop*Bwop),  Czop',     (Pop*Aop)'*Top+Top'*(Pop*Aop)];
% eps on the w and z blocks only; the legacy concatenation merged all R
% parts into one R block ordered w, z, then the ODE state.
nwz = size(Bwop,2) + size(Czop,1);
Mep = blkdiag(ep*eye(nwz),zeros(Dop.dim(1,1)-nwz),zeros(Dop.dim(2,1)));
Dop = Dop + mat2opvar(Mep,Dop.dim(:,2),PIE.vars,PIE.dom);
[prog,Deop] = poslpivar(prog,Dop.dim,dd2,opt0);
prog = lpi_eq(prog,Deop+Dop,'symmetric');
nd = numel(prog.decvartable);
evalc('sol = lpisolve(prog,sopts);');
[st,info] = read_sol(sol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st,info,nd,Km] = run_copvar(PIE,gam,dd1,dd2,sopts)
% The same LPI through the containers. The KYP operator is assembled with
% the container horzcat/vertcat, which keep the spaces w, z, R^n and L_2
% separate; the negativity slack is declared over exactly those spaces.
ep = 1e-4;
Tm  = opvar2copvar(PIE.T);      Am  = opvar2copvar(PIE.A);
Bw  = opvar2copvar(PIE.B1);     Cz  = opvar2copvar(PIE.C1);
Dzw = opvar2copvar(PIE.D11);
Iw = opvar2copvar(mat2opvar(eye(size(PIE.B1,2)),PIE.B1.dim(:,2),PIE.vars,PIE.dom));
Iz = opvar2copvar(mat2opvar(eye(size(PIE.C1,1)),PIE.C1.dim(:,1),PIE.vars,PIE.dom));
Ix = opvar2copvar(mat2opvar(eye(sum(PIE.T.dim(:,1))),PIE.T.dim(:,2),PIE.vars,PIE.dom));

prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
[sp,dm] = space_list(Tm,'out');
[prog,Pm] = poscopvar(prog,dm,sp,PIE.dom,pl2pm(dd1,sp));
Pm = Pm + ep*Ix;

PB = Pm*Bw;     PA = Pm*Am;
Km = [-(gam*Iw),   Dzw',      PB'*Tm;
       Dzw,        -(gam*Iz), Cz;
       Tm'*PB,     Cz',       PA'*Tm + Tm'*PA];
Km = Km + ep*blkdiag(Iw,Iz,0*Ix);           % eps on w and z, none on the state
[sp,dm] = space_list(Km,'out');
[prog,Nm] = poscopvar(prog,dm,sp,PIE.dom,pl2pm(dd2,sp));
% prog = lpi_eq_cdopvar(prog,Nm - (-Km),'symmetric');     % KYP = -N <= 0   % MMP, 09/25/2026 (was)
prog = lpi_eq_cdopvar(prog,Km - ((-1)*Nm),'symmetric');   % KYP = -N <= 0   % MMP, 09/25/2026
nd = numel(prog.decvartable);
evalc('sol = lpisolve(prog,sopts);');
[st,info] = read_sol(sol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sp,dm] = space_list(P,side)
% The container's spaces as poscopvar takes them: one cellstr of variable
% names per space, and the component counts.
if strcmp(side,'out'),  S = P.space_out;    dm = P.dim_out(:);
else,                   S = P.space_in;     dm = P.dim_in(:);
end
sp = cell(1,size(S,1));
for k = 1:size(S,1),    sp{k} = reshape(P.vars(S(k,:)),1,[]);    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = pl2pm(dd,sp)
% 'poslpivar' degrees to 'poscopvar' degrees, one entry per space: an R^n
% space has the identity basis, degree 0 in the integration variable; an
% L_2 space gets poslpivar's Z1, Z2, Z3 (the correspondence
% 'test_poscopvar_vs_poslpivar' verifies).
d1 = dd{1};     d2 = dd{2};     d3 = dd{3};
L2 = { struct('int',d1,'mult',0), ...
       struct('int',d2(1),'mult',d2(2),'joint',d2(3)), ...
       struct('int',d3(1),'mult',d3(2),'joint',d3(3)) };
deg = cell(1,numel(sp));
for k = 1:numel(sp)
    if isempty(sp{k}),  deg{k} = struct('int',0);
    else,               deg{k} = L2;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [st,info] = read_sol(sol)
% A CERTIFIED verdict: +1 feasible (pinf = 0, numerr = 0), -1 infeasible
% (pinf = 1, numerr = 0), 0 when the solver reports a numerical error and
% so certifies neither.
pinf = NaN;     ne = NaN;       fr = NaN;
if isfield(sol,'solinfo') && isfield(sol.solinfo,'info')
    q = sol.solinfo.info;
    if isfield(q,'pinf'),       pinf = q.pinf;      end
    if isfield(q,'numerr'),     ne = q.numerr;      end
    if isfield(q,'feasratio'),  fr = q.feasratio;   end
end
st = 0;
if ne==0 && pinf==0,    st = +1;    end
if ne==0 && pinf==1,    st = -1;    end
info = sprintf('pinf=%g numerr=%g fr=%+.2f',pinf,ne,fr);
end


function s = ternary(c,a,b)
if c,   s = a;  else,   s = b;  end
end
