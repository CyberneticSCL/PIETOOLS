function C = pielr_bench_cases(tier)                                        % CC, 09/23/2026
% PIELR_BENCH_CASES  The benchmark catalogue: what the low-rank certifier is
% exercised on.
%
% WHY A CATALOGUE AND NOT ONE BENCHMARK.  Every measurement this package has
% recorded was taken on one program -- a 2-D reaction-diffusion at one degree
% setting -- and a change measured there says nothing about a change's effect
% anywhere else.  The cases below vary the three axes that actually move the
% problem: the EXECUTIVE (a feasibility LPI vs one with an objective), the
% MODEL TYPE (1-D PDE, 2-D PDE, DDE -- which is the only family here with a
% nonzero ODE block), and the SETTINGS (basis degree).
%
% NO REFERENCE VALUE IS STORED HERE.  pielr_bench solves the identical
% assembled program with the interior-point method in the same session and
% judges it with the same gate, so the reference is measured, never quoted.
% The .note fields carry literature values as CONTEXT only; nothing scores
% against them.
%
% TIERS, because the 2-D cases cost hours and the rest cost seconds:
%   1  1-D PDE and DDE            (seconds each)
%   2  adds 2-D PDE at low degree (minutes each)
%   3  adds 2-D PDE at degree 4   (hours; the historical benchmark)
%
% INPUT   tier  1, 2 or 3 (default 1)
% OUTPUT  C     struct array with fields
%           .name .family .lpi .settings .tier .load (fn handle -> PIE) .note

if nargin<1 || isempty(tier), tier = 1; end
C = struct('name',{},'family',{},'lpi',{},'settings',{},'tier',{}, ...
           'load',{},'note',{});

% ===================== tier 1: 1-D PDE, stability ========================
% The reaction-diffusion family is the one the package's measured "rank 2"
% claims come from; included at two operating points because the rank is
% expected to climb as the operating point approaches the certifiable
% boundary, and a benchmark that only probes the easy point cannot see that.
C = add(C,'rd1d-lam0.5','1d-pde','stability','light',1, ...
    @() rd_pie(1,0.50*pi^2), 'x_t = x_ss + lam x, Dirichlet; stable for lam < pi^2');
C = add(C,'rd1d-lam0.9','1d-pde','stability','light',1, ...
    @() rd_pie(1,0.90*pi^2), 'same family, near the boundary');
C = add(C,'rd1d-n2','1d-pde','stability','light',1, ...
    @() rd_pie(2,0.50*pi^2), 'two replicated states: the 2n rank law should bite here');
C = add(C,'advdiff1d','1d-pde','stability','light',1, ...
    @() advdiff_pie(3,2), 'x_t = x_ss + 3x_s + 2x: NOT self-adjoint, so primal and dual differ');
C = add(C,'advdiff1d-dual','1d-pde','stability-dual','light',1, ...
    @() advdiff_pie(3,2), 'the dual arm on the same non-self-adjoint system');

% library stability examples: structurally different 1-D systems, several
% with an ODE component, which the reaction-diffusion family does not have
C = add(C,'lib-transport','1d-pde','stability','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Transport_Eq'), 'hyperbolic, not parabolic');
C = add(C,'lib-heat-ode-das','1d-pde','stability','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das'), 'PDE coupled to an ODE: nonzero n0 block');
C = add(C,'lib-beam-eb','1d-pde','stability','heavy',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq'), 'fourth order; the Example-21 family where max|Dop| collapsed');
C = add(C,'lib-wave-damped','1d-pde','stability','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Damped'), 'second order in time');

% ===================== tier 1: 1-D PDE, l2gain ===========================
% These carry an objective, so they score on gamma against the interior-point
% gamma rather than on a binary gate -- the only cases in the suite that give
% a continuous error measure.
C = add(C,'gain-transport','1d-pde','l2gain','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance'), 'library note: gamma ~ 0.5162 at light');
C = add(C,'gain-transport-dual','1d-pde','l2gain-dual','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance'), 'dual KYP on the same system');
C = add(C,'gain-heat-dist','1d-pde','l2gain','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance'), 'library note: gamma ~ 0.2083 at veryheavy');
C = add(C,'gain-reacdiff','1d-pde','l2gain','light',1, ...
    @() lib_pie('PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance'), 'library note: gamma ~ 8.21 near lam = pi^2');

% ===================== tier 1: DDE =======================================
% The only family here whose PIE has a substantial ODE block.  The 2-D
% benchmark prints dims [n0 nx ny n2] = [0 0 0 1], so nothing measured so far
% touches an operator with a real finite-dimensional part, and the square-root
% mechanism behind the rank claims says nothing about one.  Each is tiny.
C = add(C,'dde-scalar','dde','stability','light',1, ...
    @() dde_pie(struct('A0',0,'Ai',{{-1}},'tau',1.4)), 'xdot = -x(t-tau); literature taumax = 1.5707');
C = add(C,'dde-2state','dde','stability','light',1, ...
    @() dde_pie(struct('A0',[0 1;-2 .1],'Ai',{{[0 0;1 0]}},'tau',1.5)), 'literature taumax = 1.71785');
C = add(C,'dde-2delay','dde','stability','light',1, ...
    @() dde_pie(struct('A0',-2,'Ai',{{2.99,-1}},'tau',[1 2])), 'two delays; literature bmax = 3');
C = add(C,'dde-2state-dual','dde','stability-dual','light',1, ...
    @() dde_pie(struct('A0',[-2 0;0 -.9],'Ai',{{[-1 0;-1 -1]}},'tau',6)), 'literature stable for tau < 6.17258');

% ===================== tier 2: 2-D PDE, affordable =======================
% Operating points chosen from an IPM-only probe (interior-point only, ~120 s
% per point at degree 3) rather than guessed, because a low-rank run there
% costs 8.7x the solve it is being compared against:
%   deg 3, frac 0.10 : IPM certifies, rel 9.459e-07   <- the usable case
%   deg 3, frac 0.25 : IPM fails,     rel 5.002e-06
%   deg 3, frac 0.50 : IPM fails,     rel 7.094e-05
%   deg 2, frac 0.50 : IPM fails,     rel 4.935e+00
% So degree 2 and degree 3 at frac 0.50 are NEGATIVE CONTROLS -- the basis is
% too coarse for the operating point and neither solver can certify.  They
% earn their place: a suite with no case that must fail cannot tell a method
% that is working from a gate that accepts anything.
C = add(C,'rd2d-deg3-f010','2d-pde','stability',[],2, ...
    @() rd2d_pie(1,0.10*2*pi^2), 'degree 3 at a tenth of the boundary; IPM certifies at rel 9.459e-07');
C = add(C,'rd2d-deg3-f010-dual','2d-pde','stability-dual',[],2, ...
    @() rd2d_pie(1,0.10*2*pi^2), 'the dual arm on the same 2-D program');
C = add(C,'rd2d-deg2','2d-pde','stability',[],2, ...
    @() rd2d_pie(1,0.50*2*pi^2), 'NEGATIVE CONTROL: basis too coarse, IPM fails at rel 4.935');
C = add(C,'rd2d-deg3','2d-pde','stability',[],2, ...
    @() rd2d_pie(1,0.50*2*pi^2), 'NEGATIVE CONTROL: marginal, IPM fails at rel 7.094e-05');

% ===================== tier 3: the historical benchmark ==================
C = add(C,'rd2d-deg4','2d-pde','stability',[],3, ...
    @() rd2d_pie(1,0.50*2*pi^2), 'the program every existing measurement was taken on');

% degree settings for the 2-D cases, which take a struct rather than a name
for k = 1:numel(C)
    if strcmp(C(k).family,'2d-pde') && isempty(C(k).settings)
        d = 4;
        if contains(C(k).name,'deg2'), d = 2; elseif contains(C(k).name,'deg3'), d = 3; end
        C(k).settings = set2d_deg(d,[]);
    end
end
C = C([C.tier] <= tier);
end

% =========================================================================
function C = add(C,name,family,lpi,settings,tier,load,note)
C(end+1) = struct('name',name,'family',family,'lpi',lpi,'settings',settings, ...
                  'tier',tier,'load',load,'note',note);
end

% ---- loaders -------------------------------------------------------------
function PIE = rd_pie(n,lam)
pvar s t
phi = pde_var('state',n,s,[0,1]);
sys = [diff(phi,t,1)==diff(phi,s,2)+lam*phi;
       subs(phi,s,0)==zeros(n,1);
       subs(phi,s,1)==zeros(n,1)];
PIE = convert(sys);
end

function PIE = advdiff_pie(c,lam)
pvar s t
phi = pde_var('state',1,s,[0,1]);
sys = [diff(phi,t,1)==diff(phi,s,2)+c*diff(phi,s,1)+lam*phi;
       subs(phi,s,0)==0;  subs(phi,s,1)==0];
PIE = convert(sys);
end

function PIE = lib_pie(fname)
% The example FILES are called directly: examples_PDE_library_PIETOOLS ends in
% an unguarded input() prompt (line 976) and cannot be driven non-interactively.
f = str2func(fname);
[PDE,~] = f(0,{});
PDE = initialize_PIETOOLS_PDE(PDE,true);
PIE = convert(PDE,'pie');
end

function PIE = dde_pie(D)
% The DDE library is a script of commented-out blocks ("simply uncomment the
% variable definitions"), so it cannot be indexed; the structs are written out
% in the catalogue instead.
DDE = D;
DDE = initialize_PIETOOLS_DDE(DDE);
% one output only when an output type is named: with two, convert_PIETOOLS_DDE
% errors "At most one output is returned when an output system type is specified"
PIE = convert_PIETOOLS_DDE(DDE,'pie');
end

function PIE = rd2d_pie(n,lam)
PIE = nb_rd2d(n,lam);
end
