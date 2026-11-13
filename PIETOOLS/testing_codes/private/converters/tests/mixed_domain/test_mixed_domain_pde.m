function R = test_mixed_domain_pde(parts,names)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = TEST_MIXED_DOMAIN_PDE(PARTS,NAMES) tests PIETOOLS on coupled systems
% whose states live on spaces of different dimension (R, L2[s1], L2[s2],
% L2[s1,s2], L2[s1,s2,s3]), the systems of MIXED_DOMAIN_SYSTEMS, in four
% layers:
%   'convert'  the parser and converter: for each supported system, T, A,
%              Cz, K and pde2pie of the PIE against the PDE written by hand
%              (T xf = x, A xf = right-hand side, z, K xf = 0), on states
%              satisfying the BCs and, without auxiliary states, on random
%              fundamental states (BCs of T xf, D^alpha T xf = xf); each
%              operator applied by its kernels (apply_sopvar), exactly on
%              polynomials. A system that must be rejected passes if the
%              conversion errors with its keyword;
%   'ops'      the 'copvar' operations on these PIE operators: A*T, T'*A,
%              T+A, T'*A+A'*T against sequential application, and T', A'
%              against the inner product (integrated exactly);
%   'exact'    the exact baselines of 'exec': stability margins, H-infinity
%              gains and H2 norms from separation of variables, each
%              checked against the collocation model MIXED_DOMAIN_COLLOC;
%   'exec'     PIETOOLS_stability_sop ('P') at 0.9 and 1.1 times each exact
%              margin (certify / must not certify), PIETOOLS_Hinf_gain_sop
%              ('Q', 'P') against the exact gains, and PIETOOLS_H2_norm_sop ('o',
%              'oco') against the exact H2 norms; lpisettings_sop('light'),
%              MOSEK when installed. About 1 h;
%   'piesim'   PIESIM (BDF4, N = 16) against the exact solutions of the
%              PIESIM 2-D library examples 8-11, the collocation model, and
%              dde23 on library PDE 43.
% Prints one line per check.
%
% INPUT
% - parts:  (optional) char or cellstr, the layers to run (default
%           {'convert','ops','exact'}, about 3 min; 'all' runs every layer);
% - names:  (optional) cellstr, the systems of 'convert' and 'ops' (default
%           every system of the catalogue).
%
% OUTPUT
% - R:      struct array, one element per check, with fields part, name,
%           status ('pass', 'FAIL', 'KNOWN' for a recorded limitation,
%           'rejected' for an expected rejection) and detail (char).
%
% NOTES
% The folder must be on the path: testing_codes/private/ is off the
% default path, e.g.
%   addpath(fullfile(<PIETOOLS>,'testing_codes','private','converters','tests','mixed_domain'))
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - test_mixed_domain_pde
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 11/12/2026

% Layers and systems to run
if nargin<1 || isempty(parts), parts = {'convert','ops','exact'}; end
parts = cellstr(parts);
if any(strcmp(parts,'all')), parts = {'convert','ops','exact','exec','piesim'}; end
catl = mixed_domain_systems('list');
if nargin<2 || isempty(names), names = {catl.name}; end
warning('off','sopvar:noncanonicalMultiplier');  warning('off','sdopvar:noncanonicalMultiplier');
R = struct('part',{},'name',{},'status',{},'detail',{});
for k=1:numel(parts)
    switch parts{k}
        case 'convert', R = [R, run_convert(names)];                       %#ok<AGROW>
        case 'ops',     R = [R, run_ops(names)];                           %#ok<AGROW>
        case 'exact',   R = [R, run_exact()];                              %#ok<AGROW>
        case 'exec',    R = [R, run_exec(names)];                          %#ok<AGROW>
        case 'piesim',  R = [R, run_piesim()];                             %#ok<AGROW>
        otherwise, error('test_mixed_domain_pde:part','Unknown part %s.',parts{k});
    end
end
% Summary
st = {R.status};
fprintf('\n%d checks: %d pass, %d rejected as required, %d KNOWN, %d FAIL.\n',numel(R), ...
    nnz(strcmp(st,'pass')),nnz(strcmp(st,'rejected')),nnz(strcmp(st,'KNOWN')),nnz(strncmp(st,'FAIL',4)));
end



%% ------------------------------------------------------------------------
%% Layer 1: parser and converter
%% ------------------------------------------------------------------------
function R = run_convert(names)
% R = RUN_CONVERT(NAMES) converts each system and checks its PIE against
% the hand-written semantics (see the help of the file).
%
% INPUT
% - names:  cellstr, the systems.
%
% OUTPUT
% - R:      struct array of results (part 'convert').
tol = 1e-10;
known = known_limitations();
R = struct('part',{},'name',{},'status',{},'detail',{});
for i=1:numel(names)
    S = mixed_domain_systems(names{i});
    r = struct('part','convert','name',names{i},'status','','detail','');
    % Convert, and classify a conversion error ('silent', not evalc: evalc
    % in a case loop crashes the R2025b interpreter, access violation)
    try
        PIE = convert(S.PDE,{'pie','silent'});
    catch ME
        msg = ME.message;
        if ~isempty(S.expect) && contains(lower(msg),S.expect)
            r.status = 'rejected';
        elseif isfield(known,names{i}) && contains(msg,known.(names{i}))
            r.status = 'KNOWN';
        else
            r.status = 'FAIL';
        end
        r.detail = msg;
        R(end+1) = report(r);                                              %#ok<AGROW>
        continue
    end
    if ~isempty(S.expect)
        r.status = 'FAIL';  r.detail = 'converted, but must be rejected';
        R(end+1) = report(r);                                              %#ok<AGROW>
        continue
    end
    % Semantic checks
    try
        e = check_pie(S,PIE,2);
        v = [e.eT e.eA e.eZ e.eK e.eC e.eBC e.eD e.eAB];
        r.status = ternary(all(v(~isnan(v))<tol),'pass','FAIL');
        r.detail = sprintf('blocks %d (aux %d) | T %.1e A %.1e Cz %.1e K %.1e C %.1e | BC %.1e D %.1e A(rand) %.1e', ...
            e.nblk,e.naux,e.eT,e.eA,e.eZ,e.eK,e.eC,e.eBC,e.eD,e.eAB);
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
end

function K = known_limitations()
% K = KNOWN_LIMITATIONS() returns, per system, the text of the error that
% records a measured limitation of the converter (11/12/2026).
K = struct('S19_2d_integral_in_1d_bc','cannot be separated into conditions along a single spatial direction');
end

function e = check_pie(S,PIE,nrep)
% E = CHECK_PIE(S,PIE,NREP) checks the PIE of system S, NREP random test
% functions per check.
%
% INPUT
% - S:      'struct' of MIXED_DOMAIN_SYSTEMS;
% - PIE:    'pie_struct', the converted PIE;
% - nrep:   double, the number of test functions.
%
% OUTPUT
% - e:      'struct' with the relative errors eT, eA, eZ, eK, eC (states
%           built to satisfy the BCs) and eBC, eD, eAB (random fundamental
%           states, NaN with auxiliary states), and nblk, naux.
T = PIE.T;  A = PIE.A;
vreg = T.vars;  dom = T.dom;
pvars = PIE.space.vars;                    % column order of pde2pie.x.diff
idx = PIE.pde2pie.x.idx(:,1);  dmat = PIE.pde2pie.x.diff;
nb = numel(idx);
% Variables of each block; an auxiliary block lacks some of its state's
% variables and holds the integral of the state over them
bvars = cell(nb,1);  drop = cell(nb,1);
for b=1:nb
    bvars{b} = vreg(T.space_in(b,:));
    drop{b} = setdiff(S.svars{idx(b)},bvars{b});
end
e = struct('nblk',nb,'naux',nnz(~cellfun(@isempty,drop)),'eT',NaN,'eA',NaN,'eZ',NaN,'eK',NaN, ...
           'eC',NaN,'eBC',NaN,'eD',NaN,'eAB',NaN);
hasK = nonempty(PIE.K);
% Fundamental inputs: each must be a PDE input itself
widx = [];
if nonempty(PIE.Bw) || nonempty(PIE.Tw) || nonempty(PIE.Dzw)
    widx = PIE.pde2pie.w.idx(:);
    if any(PIE.pde2pie.w.diff(:)) || any(PIE.pde2pie.w.loc(:))
        error('A fundamental input is a derivative or boundary value of a PDE input (not handled).');
    end
end
nz = numel(S.z(repmat({polynomial(0)},1,numel(S.svars)),repmat({polynomial(0)},1,numel(S.wvars))));
% (A) states built by hand to satisfy the BCs
eT = 0; eA = 0; eZ = 0; eK = 0; eC = 0;
for rep=1:nrep
    rng(100+rep);
    c = randn(1,40);
    W = cell(1,numel(S.wvars));
    for j=1:numel(W), W{j} = randpoly(S.wvars{j},2,1); end
    X = S.xs(c,W);  F = S.rhs(X,W);  bres = S.bc(X,W);
    for j=1:numel(bres)
        if pmax(bres{j})>1e-12*max(1,maxc(X))
            error('The hand-built state violates BC %d (test error).',j);
        end
    end
    % xf_b = F_b D^{diff_b} x_k; rows of T xf: x_k or F_b x_k; of A xf: the
    % right-hand side, or F_b times it
    xf = cell(nb,1);  Trow = cell(nb,1);  Arow = cell(nb,1);
    for b=1:nb
        k = idx(b);
        xf{b} = Fop(Dop(X{k},pvars,dmat(b,:)),drop{b},vreg,dom);
        Trow{b} = Fop(X{k},drop{b},vreg,dom);
        Arow{b} = Fop(F{k},drop{b},vreg,dom);
    end
    wf = W(widx);
    Txf = addc(apply_cop(T,xf),apply_cop(PIE.Tw,wf,nb));
    Axf = addc(apply_cop(A,xf),apply_cop(PIE.Bw,wf,nb));
    for b=1:nb
        eT = max(eT,rel(Txf{b},Trow{b}));
        eA = max(eA,rel(Axf{b},Arow{b}));
    end
    if nz>0
        Zr = S.z(X,W);
        Zx = addc(apply_cop(PIE.Cz,xf,nz),apply_cop(PIE.Dzw,wf,nz));
        for j=1:nz, eZ = max(eZ,rel(Zx{j},Zr{j})); end
    end
    if hasK
        Kxf = apply_cop(PIE.K,xf);
        for j=1:numel(Kxf), eK = max(eK,pmax(Kxf{j})/max(1,maxc(xf))); end
    end
    % pde2pie.x.C: one column per block (its state's space), diagonal
    Cm = PIE.pde2pie.x.C;
    for b=1:nb
        v = cell(size(Cm.C,2),1);
        for kk=1:size(Cm.C,2), v{kk} = polynomial(zeros(Cm.dim_in(kk),1)); end
        v{b} = Dop(X{idx(b)},pvars,dmat(b,:));
        y = apply_cop(Cm,v);
        eC = max(eC,rel(y{b},xf{b}));
    end
end
e.eT = eT;  e.eA = eA;  e.eK = eK;  e.eC = eC;
if nz>0, e.eZ = eZ; end
if hasK, return, end
% (B) random fundamental states and inputs, no auxiliary blocks
eBC = 0; eD = 0; eAB = 0;
for rep=1:nrep
    rng(200+rep);
    xf = cell(nb,1);
    for b=1:nb, xf{b} = randpoly(bvars{b},3,T.dim_in(b)); end
    W = cell(1,numel(S.wvars));
    for j=1:numel(W), W{j} = randpoly(S.wvars{j},2,1); end
    wf = W(widx);
    x = addc(apply_cop(T,xf),apply_cop(PIE.Tw,wf,nb));
    Axf = addc(apply_cop(A,xf),apply_cop(PIE.Bw,wf,nb));
    X = cell(1,numel(S.svars));
    for b=1:nb, X{idx(b)} = x{b}; end
    bres = S.bc(X,W);
    for j=1:numel(bres), eBC = max(eBC,pmax(bres{j})/max(1,maxc(x))); end
    F = S.rhs(X,W);
    for b=1:nb
        k = idx(b);
        eD = max(eD,rel(Dop(X{k},pvars,dmat(b,:)),xf{b}));
        eAB = max(eAB,rel(Axf{b},F{k}));
    end
end
e.eBC = eBC;  e.eD = eD;  e.eAB = eAB;
end



%% ------------------------------------------------------------------------
%% Layer 2: copvar operations
%% ------------------------------------------------------------------------
function R = run_ops(names)
% R = RUN_OPS(NAMES) checks the 'copvar' operations on the operators T, A
% of each supported system (see the help of the file).
%
% INPUT
% - names:  cellstr, the systems (those to be rejected are skipped).
%
% OUTPUT
% - R:      struct array of results (part 'ops').
tol = 1e-10;
R = struct('part',{},'name',{},'status',{},'detail',{});
for i=1:numel(names)
    S = mixed_domain_systems(names{i});
    if ~isempty(S.expect) || isfield(known_limitations(),names{i}), continue, end
    r = struct('part','ops','name',names{i},'status','','detail','');
    try
        PIE = convert(S.PDE,{'pie','silent'});
        T = PIE.T;  A = PIE.A;
        AT = A*T;  Tt = T';  At = A';  TtA = Tt*A;  TA = T+A;  Lop = Tt*A+At*T;
        eM = 0; eMt = 0; eAdj = 0; eS = 0; eL = 0;
        for rep=1:2
            rng(300+rep);
            f = rand_blocks(T,'in');  g = rand_blocks(T,'out');
            Tf = apply_cop(T,f);  Af = apply_cop(A,f);  TtAf = apply_cop(Tt,Af);
            eM  = max(eM, cmpc(apply_cop(AT,f),apply_cop(A,Tf)));
            eMt = max(eMt,cmpc(apply_cop(TtA,f),TtAf));
            eS  = max(eS, cmpc(apply_cop(TA,f),addc(Tf,Af)));
            eL  = max(eL, cmpc(apply_cop(Lop,f),addc(TtAf,apply_cop(At,Tf))));
            i1 = ip(Tf,g,T,'out');  i2 = ip(Af,g,A,'out');
            e1 = i1-ip(f,apply_cop(Tt,g),T,'in');
            e2 = i2-ip(f,apply_cop(At,g),A,'in');
            eAdj = max(eAdj,max(abs([e1 e2]))/max([1,abs(i1),abs(i2)]));
        end
        r.status = ternary(all([eM eMt eAdj eS eL]<tol),'pass','FAIL');
        r.detail = sprintf('A*T %.1e  T''*A %.1e  adjoint %.1e  T+A %.1e  T''A+A''T %.1e',eM,eMt,eAdj,eS,eL);
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
end



%% ------------------------------------------------------------------------
%% Exact baselines
%% ------------------------------------------------------------------------
function [R,E] = run_exact()
% [R,E] = RUN_EXACT() computes the exact baselines of 'exec' and checks
% each against the collocation model.
%
% OUTPUT
% - R:      struct array of results (part 'exact');
% - E:      'struct' of the baselines (see EXACT_BASELINES).
E = exact_baselines();
R = struct('part',{},'name',{},'status',{},'detail',{});
% Collocation self-test
st = mixed_domain_colloc('selftest',[],24);
r = struct('part','exact','name','colloc selftest','status', ...
    ternary(st.asc && st.eD<1e-12 && st.eD2<1e-9 && st.eW<1e-14,'pass','FAIL'), ...
    'detail',sprintf('d/dx %.1e, d2/dx2 %.1e, quadrature %.1e',st.eD,st.eD2,st.eW));
R(end+1) = report(r);
% Rightmost eigenvalue of the collocation model at each exact margin is 0,
% and at the S9 point a = 0.5 it is the root of the characteristic equation
for i=1:numel(E.margins)
    m = E.margins(i);
    M = mixed_domain_colloc(m.sys,m.atstar,24);
    lam = real(M.eig(1));
    r = struct('part','exact','name',[m.sys,' margin'],'status',ternary(abs(lam-m.lam)<1e-7,'pass','FAIL'), ...
        'detail',sprintf('%s* = %.10g: collocation max Re %.2e, exact %.2e',m.par,m.star,lam,m.lam));
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
% H-infinity: the scalar-input gains against the collocation model
for k={'h1','h3'}
    h = E.(k{1});
    r = struct('part','exact','name',[k{1},' ',h.sys],'status',ternary(abs(h.colloc/h.gam-1)<1e-8,'pass','FAIL'), ...
        'detail',sprintf('gain %.10g at omega %.3g; collocation %.10g',h.gam,h.omega,h.colloc));
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
% H2: the quadrature value against the collocation Gramian
for k={'n1','n3'}
    h = E.(k{1});
    r = struct('part','exact','name',['H2 ',h.sys],'status',ternary(abs(h.colloc/h.h2-1)<1e-8,'pass','FAIL'), ...
        'detail',sprintf('H2 norm %.10g by quadrature; collocation %.10g',h.h2,h.colloc));
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
% Control C1: the gain against the double series (1e-9) and the
% collocation; the H2 norm against the collocation. The uniform source does
% not vanish on the boundary, so the solution has corner singularities and
% the collocation converges algebraically: measured 1.1e-9 (gain) and
% 9.6e-9 (H2) at N = 24 on 11/12/2026, against 1e-10 for the smooth cases;
% tolerance 1e-7, not 1e-8.
h = E.h5;  n = E.n5;
r = struct('part','exact','name',['h5 ',h.sys],'status', ...
    ternary(abs(h.series/h.gam-1)<1e-9 && abs(h.colloc/h.gam-1)<1e-7,'pass','FAIL'), ...
    'detail',sprintf('gain %.10g; double series %.10g; collocation %.10g',h.gam,h.series,h.colloc));
R(end+1) = report(r);
r = struct('part','exact','name',['H2 ',n.sys],'status',ternary(abs(n.colloc/n.h2-1)<1e-7,'pass','FAIL'), ...
    'detail',sprintf('H2 norm %.10g by quadrature; collocation %.10g',n.h2,n.colloc));
R(end+1) = report(r);
end

function E = exact_baselines()
% E = EXACT_BASELINES() returns the exact stability margins and H-infinity
% gains of the executive tests, by separation of variables.
%
% OUTPUT
% - E:      'struct' with fields
%             margins  struct array: sys, par (the parameter), star (its
%                      critical value), atstar ('struct' of parameters at
%                      the margin), lam (rightmost eigenvalue there, 0 at a
%                      margin) and base ('struct' of the other parameters);
%             h1..h5   struct: sys, par, gam (the gain), omega (where it is
%                      attained), colloc (h1, h3, h5: collocation value),
%                      series (h5: the gain by a double series);
%             n1,n3,n5 struct: sys, par, h2 (the H2 norm), colloc;
%             tau43    the delay margin of library PDE 43 (sign -1).
% h5 and n5 belong to the control C1 (one 2-D heat state, no lower-
% dimensional or auxiliary state).
% S4 Wentzell, c=1: lam = r - m^2 pi^2 - k^2, k tan k = c
k0 = fzero(@(k) k*tan(k)-1,[0.1 1.5]);
% S9, c=1: lam + c mu coth mu - a = 0, mu^2 = lam; a* = c; at a = 0.5 the
% rightmost root is -k^2 with -k^2 + k cot k - 0.5 = 0
k9 = fzero(@(k) -k^2+k*cot(k)-0.5,[0.05 1.5]);
mg = @(sys,par,star,base,lam) struct('sys',sys,'par',par,'star',star, ...
        'atstar',setfield(base,par,star),'lam',lam,'base',base);
E.margins = [mg('S3_edge1d_heat2d','r',pi^2,struct(),0);
             mg('S4_wentzell','r',pi^2+k0^2,struct('c',1),0);
             mg('S6_indomain_1d_2d','g',pi^2,struct(),0);
             mg('S7_ode_2d_indomain','g',1,struct('a',1,'b',1),0);
             mg('S8_two_edges','r3',pi^2,struct(),0);
             mg('S9_ode_2d_flux','a',1,struct('c',1),0);
             mg('S9_ode_2d_flux','a',0.5,struct('c',1),-k9^2)];
E.margins(end).star = 0.5;                  % a check point, not a margin
% Library PDE 43 with the sign of its reference: lam = 0.9 - exp(-lam tau)
w43 = sqrt(1-0.81);  E.tau43 = atan2(w43,0.9)/w43;
% H-infinity gains: sup over the imaginary axis
wg = [0, logspace(-3,3,3000)];
E.h1 = hinf_scalar(@(s) 1./(s-0.5+mucoth(s)),wg);
E.h1.sys = 'S9_ode_2d_flux';  E.h1.par = struct('a',0.5,'c',1);
E.h3 = hinf_scalar(@(s) 1./(s+1).*hint(s),wg);
E.h3.sys = 'S24_hinf_ode_2d_edge';  E.h3.par = struct('a',1);
E.h2 = hinf_norm(@(s) series_indomain(s,1),wg);
E.h2.sys = 'S23_hinf_indomain';  E.h2.par = struct('g',1);
E.h4 = hinf_norm(@(s) series_wentzell(s,1,0),wg);
E.h4.sys = 'S25_hinf_wentzell';  E.h4.par = struct('c',1,'r',0);
% Scalar-input gains from the collocation model, G = C (sI - L)^{-1} B
wc = [0, logspace(-3,3,400)];
M = mixed_domain_colloc('S9_ode_2d_flux',struct('a',0.5,'c',1),24);
bf = zeros(numel(M.idxI)+numel(M.idxB),1);  bf(M.state(1).idx) = 1;
E.h1.colloc = colloc_gain(M,bf,bf',wc);
% H2 norm of the same input-output map: (1/pi) int_0^inf |G(i w)|^2 dw,
% against C X C' with L X + X L' + B B' = 0 on the collocation model
E.n1 = struct('sys','S9_ode_2d_flux','par',struct('a',0.5,'c',1), ...
              'h2',h2_quad(@(s) 1./(s-0.5+mucoth(s))),'colloc',colloc_h2(M,bf,bf'));
M = mixed_domain_colloc('S2_ode_heat2d_edge',struct('a',1),24);
bf = zeros(numel(M.idxI)+numel(M.idxB),1);  bf(M.state(1).idx) = 1;
cf = zeros(1,numel(bf));  cf(M.state(2).idx) = kron(M.w(:)',M.w(:)');
E.h3.colloc = colloc_gain(M,bf,cf,wc);
E.n3 = struct('sys','S24_hinf_ode_2d_edge','par',struct('a',1), ...
              'h2',h2_quad(@(s) hint(s)./(s+1)),'colloc',colloc_h2(M,bf,cf));
% Control C1: modes 2 sin(m pi s1) sin(n pi s2), m, n odd, carry the
% residue a_m a_n, a_m = 8/(m^2 pi^2) > 0, so the impulse response is f(t)^2
% with f(t) = sum_{m odd} a_m exp(-m^2 pi^2 t). All residues are positive,
% so the gain is G(0) = int f^2 dt, and the H2 norm is sqrt(int f^4 dt).
m = (1:2:40001)';  am = 8./(m.^2*pi^2);  lm = m.^2*pi^2;
f = @(t) sum(am.*exp(-lm.*t),1);
E.h5 = struct('sys','C1_heat2d_scalar_w','par',struct(),'omega',0, ...
              'gam',integral(@(t) f(t).^2,0,Inf,'AbsTol',1e-15,'RelTol',1e-12));
mm = (1:2:2001)';  A5 = 8./(mm.^2*pi^2);                % tail below 1e-10
E.h5.series = sum(sum((A5*A5')./((mm.^2+mm'.^2)*pi^2)));
E.n5 = struct('sys','C1_heat2d_scalar_w','par',struct(), ...
              'h2',sqrt(integral(@(t) f(t).^4,0,Inf,'AbsTol',1e-15,'RelTol',1e-12)));
M = mixed_domain_colloc('C1_heat2d_scalar_w',struct(),24);
bf = zeros(numel(M.idxI)+numel(M.idxB),1);  bf(M.state(1).idx) = 1;
cf = zeros(1,numel(bf));  cf(M.state(1).idx) = kron(M.w(:)',M.w(:)');
E.h5.colloc = colloc_gain(M,bf,cf,wc);
E.n5.colloc = colloc_h2(M,bf,cf);
end

function v = h2_quad(G)
% V = H2_QUAD(G) returns the H2 norm sqrt((1/pi) int_0^inf |G(i w)|^2 dw) of
% the scalar transfer function G (vectorized), by adaptive quadrature.
v = sqrt(integral(@(w) abs(G(1i*w)).^2,0,Inf,'AbsTol',1e-14,'RelTol',1e-12)/pi);
end

function v = colloc_h2(M,bf,cf)
% V = COLLOC_H2(M,BF,CF) returns the H2 norm of the collocation model M
% from the input column BF to the output row CF: sqrt(C X C') with
% L X + X L' + B B' = 0.
B = bf(M.idxI);  C = cf(M.idxI)+cf(M.idxB)*M.G;
X = sylvester(M.L,M.L',-B*B');
v = sqrt(real(C*X*C'));
end

function v = mucoth(s)
% V = MUCOTH(S) returns sqrt(s) coth(sqrt(s)), 1 at s = 0.
mu = sqrt(s);  v = ones(size(s));  nz = abs(mu)>1e-8;
v(nz) = mu(nz).*coth(mu(nz));
end

function v = hint(s)
% V = HINT(S) returns int_0^1 sinh(mu(1-x))/sinh(mu) dx = (cosh mu - 1)/
% (mu sinh mu) = tanh(mu/2)/mu, mu = sqrt(s), 1/2 at s = 0; the tanh form
% does not overflow at large |s|.
mu = sqrt(s);  v = 0.5*ones(size(s));  nz = abs(mu)>1e-6;
v(nz) = tanh(mu(nz)/2)./mu(nz);
end

function H = hinf_scalar(G,wg)
% H = HINF_SCALAR(G,WG) returns sup |G(i omega)| over the grid WG,
% refined by fminbnd, as 'struct' with fields gam and omega.
H = hinf_norm(@(s) abs(G(s)),wg);
end

function H = hinf_norm(nrm,wg)
% H = HINF_NORM(NRM,WG) returns sup NRM(i omega) over the grid WG, refined
% by fminbnd between the neighbours of the largest grid value.
g = arrayfun(@(w) nrm(1i*w),wg);
[gm,k] = max(g);  om = wg(k);
if k>1 && k<numel(wg)
    [wo,fv] = fminbnd(@(w) -nrm(1i*w),wg(k-1),wg(k+1),optimset('TolX',1e-12));
    if -fv>gm, gm = -fv; om = wo; end
end
H = struct('gam',gm,'omega',om);
end

function v = series_indomain(s,g)
% V = SERIES_INDOMAIN(S,G) returns ||G(s)|| from L2[0,1]^2 to R for S23:
% mode m (orthonormal sqrt2 sin(m pi s1), constant in s2) has
% a_m = g/((s+m^2 pi^2)^2-g^2) w_m, z = sum_m int sin(m pi s1) a_m.
m = (1:2:4001)';  k = (m*pi).^2;  cm = 2./(m*pi);
v = sqrt(sum(abs(sqrt(2)*cm.*g./((s+k).^2-g^2)).^2));
end

function v = series_wentzell(s,c,r)
% V = SERIES_WENTZELL(S,C,R) returns ||G(s)|| from L2[0,1] to R for S25:
% a_m = w_m/(mu^2 + c mu coth mu), mu^2 = s + m^2 pi^2 - r.
m = (1:2:4001)';  mu2 = s+(m*pi).^2-r;  mu = sqrt(mu2);  cm = 2./(m*pi);
v = sqrt(sum(abs(sqrt(2)*cm./(mu2+c*mu.*coth(mu))).^2));
end

function g = colloc_gain(M,bf,cf,wg)
% G = COLLOC_GAIN(M,BF,CF,WG) returns max over WG of |C (iw - L)^{-1} B|
% for the collocation model M, input column BF and output row CF on the
% full node vector.
B = bf(M.idxI);  C = cf(M.idxI)+cf(M.idxB)*M.G;  n = size(M.L,1);
g = max(arrayfun(@(w) abs(C*((1i*w*eye(n)-M.L)\B)),wg));
end



%% ------------------------------------------------------------------------
%% Layer 3: executives
%% ------------------------------------------------------------------------
function R = run_exec(names)
% R = RUN_EXEC(NAMES) runs PIETOOLS_stability_sop ('P') on each system at
% 0.9 and 1.1 times its exact margin, PIETOOLS_Hinf_gain_sop ('Q', 'P') on
% the H-infinity systems, and PIETOOLS_H2_norm_sop ('o', 'oco') on three
% scalar-input systems, at lpisettings_sop('light'). A stable PIE must be
% certified, an unstable one must not; a certified gain or H2 norm must lie
% in [g*, 1.05 g*] (below the exact value g* is unsound).
%
% INPUT
% - names:  (optional) cellstr, run only the cases of these systems.
%
% OUTPUT
% - R:      struct array of results (part 'exec'); detail holds the status,
%           the certifying degree raise dw, the rows m, the decision
%           variables and the time.
E = exact_baselines();
% {system, parameters, 'stab' or 'hinf', expected: true/false or the gain}
C = {};
for f=[0.9 1.1]
    C(end+1,:) = {'S3_edge1d_heat2d',struct('r',f*pi^2),'stab',f<1};          %#ok<AGROW>
    C(end+1,:) = {'S4_wentzell',struct('c',1,'r',f*E.margins(2).star),'stab',f<1}; %#ok<AGROW>
    C(end+1,:) = {'S6_indomain_1d_2d',struct('g',f*pi^2),'stab',f<1};         %#ok<AGROW>
    C(end+1,:) = {'S7_ode_2d_indomain',struct('a',1,'b',1,'g',f),'stab',f<1};  %#ok<AGROW>
    C(end+1,:) = {'S8_two_edges',struct('r3',f*pi^2),'stab',f<1};             %#ok<AGROW>
    C(end+1,:) = {'S9_ode_2d_flux',struct('a',f,'c',1),'stab',f<1};            %#ok<AGROW>
    C(end+1,:) = {'S17_transport_edge_diffusion',struct('r2',f*pi^2),'stab',f<1}; %#ok<AGROW>
    C(end+1,:) = {'S18_lib43_interior_delay',struct('tau',f*E.tau43,'sgn',-1),'stab',f<1}; %#ok<AGROW>
    % PDE 42: mode cos(pi s1/2) of x2 obeys b' = (a - pi^2/4) b + 0.5 b(t-1),
    % stable iff a - pi^2/4 + 0.5 < 0; the ODE is driven by x2 (cascade)
    C(end+1,:) = {'S26_lib42_ode_heat_delay',struct('a',f*(pi^2/4-0.5)),'stab',f<1}; %#ok<AGROW>
end
for a=[0.5 -0.5]
    C(end+1,:) = {'S1_ode_heat1d_bc',struct('a',a),'stab',a>0};               %#ok<AGROW>
    C(end+1,:) = {'S2_ode_heat2d_edge',struct('a',a),'stab',a>0};             %#ok<AGROW>
    C(end+1,:) = {'S5_chain_ode_1d_2d',struct('a',a),'stab',a>0};             %#ok<AGROW>
end
C(end+1,:) = {'S18_lib43_interior_delay',struct('tau',0.5,'sgn',1),'stab',false};
% h5, n5: the control C1, a 2-D PIE without lower-dimensional or
% auxiliary states, separates 2-D effects from mixed-domain ones
for k={'h1','h2','h3','h4','h5'}
    h = E.(k{1});
    C(end+1,:) = {h.sys,h.par,'hinf',h.gam};                               %#ok<AGROW>
    C(end+1,:) = {h.sys,h.par,'hinfP',h.gam};                              %#ok<AGROW>
end
% H2 norm, forms 'o' and 'oco' (n1 and n3 have K, so 'c' and 'cco' are
% posed as these)
for k={'n1','n3','n5'}
    h = E.(k{1});
    C(end+1,:) = {h.sys,h.par,'h2o',h.h2};                                 %#ok<AGROW>
    C(end+1,:) = {h.sys,h.par,'h2oco',h.h2};                               %#ok<AGROW>
end
if nargin>=1 && ~isempty(names)
    C = C(ismember(C(:,1),names),:);
end
% Settings: light tier, MOSEK when installed, quiet
st = lpisettings_sop('light');
if ~isempty(which('mosekopt')), st.sos_opts.solver = 'mosek'; end
st.sop.verbose = false;
R = struct('part',{},'name',{},'status',{},'detail',{});
for i=1:size(C,1)
    r = struct('part','exec','name',sprintf('%s %s %s',C{i,1},C{i,3},jsonencode(C{i,2})),'status','','detail','');
    try
        S = mixed_domain_systems(C{i,1},C{i,2});
        PIE = initialize(convert(S.PDE,{'pie','silent'}));
        t0 = tic;
        if strcmp(C{i,3},'stab')
            [~,~,info] = PIETOOLS_stability_sop(PIE,st,'P');
            if info.certified && ~C{i,4}
                r.status = 'FAIL';  vd = 'UNSOUND: unstable PIE certified';
            elseif ~info.certified && C{i,4}
                r.status = 'FAIL';  vd = 'conservative: stable PIE not certified';
            else
                r.status = 'pass';  vd = ternary(C{i,4},'stable, certified','unstable, not certified');
            end
        else
            switch C{i,3}
                case 'hinf',    [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,'Q');
                case 'hinfP',   [~,~,gam,info] = PIETOOLS_Hinf_gain_sop(PIE,st,'P');
                case 'h2o',     [~,~,gam,~,~,info] = PIETOOLS_H2_norm_sop(PIE,st,'o');
                case 'h2oco',   [~,~,gam,~,~,info] = PIETOOLS_H2_norm_sop(PIE,st,'oco');
            end
            gap = gam/C{i,4}-1;
            if info.certified && gap<-1e-6
                r.status = 'FAIL';  vd = sprintf('UNSOUND: gam %.8g below exact %.8g',gam,C{i,4});
            elseif ~info.certified || gap>0.05
                r.status = 'FAIL';  vd = sprintf('gam %.8g, exact %.8g (gap %+.2e)',gam,C{i,4},gap);
            else
                r.status = 'pass';  vd = sprintf('gam %.8g, exact %.8g (gap %+.2e)',gam,C{i,4},gap);
            end
        end
        tt = toc(t0);
        sh = info.shape;
        r.detail = sprintf('%s | %s, dw %d, m %d, ndv %d, %.0f s',vd,info.status,info.dw,sh.m,sh.ndv,tt);
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
end



%% ------------------------------------------------------------------------
%% Layer 4: PIESIM
%% ------------------------------------------------------------------------
function R = run_piesim()
% R = RUN_PIESIM() simulates with PIESIM (PDE input, BDF order 4, dt 1e-5,
% final time 0.1) and compares the final states with
%   - the exact solutions of the PIESIM 2-D library examples 8-11 (N = 24
%     per direction), and of a 2-D heat equation with Neumann conditions
%     at both ends of s1 (one state, an auxiliary state in the PIE);
%   - the collocation model (expm in time, N = 28, interpolated to the
%     PIESIM grid) for the systems of MIXED_DOMAIN_COLLOC (N = 16);
%   - dde23 on the first mode of library PDE 43, both signs, with a
%     history twice continuously differentiable into the solution.
% A check passes if every state has relative error below TOL. PIESIM is
% called without evalc, so its progress is printed: evalc in a case loop
% crashes the R2025b interpreter (access violation).
%
% OUTPUT
% - R:      struct array of results (part 'piesim'); detail holds the
%           error of each state.
tol = 1e-6;  tf = 0.1;
mk = @(N) struct('N',N,'tf',tf,'intScheme',1,'Norder',4,'dt',1e-5,'plot','no','ploteig','no');
opts = mk(24);
R = struct('part',{},'name',{},'status',{},'detail',{});
syms st sx sy real
pvar s1 s2 t
% PIESIM 2-D library examples with exact solutions
for k=8:11
    r = struct('part','piesim','name',sprintf('PIESIM 2-D library ex %d',k),'status','','detail','');
    try
        [PDE,uinput] = examples_pde_library_PIESIM_2D(k);
        [sol,grd] = PIESIM(PDE,opts,uinput);
        err = err_exact(sol,grd,uinput.exact,state_kinds(PDE),tf);
        r.status = ternary(max(err)<tol,'pass','FAIL');
        r.detail = ['errors ',mat2str(err',3)];
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
% One 2-D state, Neumann at both ends of s1: the PIE has the auxiliary
% state int x ds1 (L2[s2]) and K; exact solution by separation
r = struct('part','piesim','name','2-D heat, Neumann-Neumann in s1','status','','detail','');
try
    clear stateNameGenerator
    x = pde_var(1,[s1;s2],[0,1;0,1]);
    PDE = initialize([diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
           subs(diff(x,s1),s1,0)==0; subs(diff(x,s1),s1,1)==0;
           subs(x,s2,0)==0; subs(x,s2,1)==0]);
    uex = cos(pi*sx)*sin(pi*sy)*exp(-2*pi^2*st)+sin(pi*sy)*exp(-pi^2*st);
    uinput = struct('ic',subs(uex,st,0),'exact',uex);
    [sol,grd] = PIESIM(PDE,opts,uinput);
    err = err_exact(sol,grd,uinput.exact,{'xy'},tf);
    r.status = ternary(max(err)<tol,'pass','FAIL');
    r.detail = ['errors ',mat2str(err',3)];
catch ME
    r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
end
R(end+1) = report(r);
% One 2-D state with an integral BC in s1 (S31): no auxiliary state; exact
% solution cos(pi s1) sin(pi s2) exp(-2 pi^2 t)
r = struct('part','piesim','name','S31_integral_bc_2d','status','','detail','');
try
    S = mixed_domain_systems('S31_integral_bc_2d');
    uex = cos(pi*sx)*sin(pi*sy)*exp(-2*pi^2*st);
    uinput = struct('ic',subs(uex,st,0),'exact',uex);
    [sol,grd] = PIESIM(initialize(S.PDE),opts,uinput);
    err = err_exact(sol,grd,uinput.exact,{'xy'},tf);
    r.status = ternary(max(err)<tol,'pass','FAIL');
    r.detail = ['errors ',mat2str(err',3)];
catch ME
    r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
end
R(end+1) = report(r);
% The discrete generator Atotal = Tcheb^{-1} Acheb of PIESIM (N = 16): its
% rightmost eigenvalue against the exact one, for ordinary BCs (Neumann at
% s1 = 0, Dirichlet elsewhere), an integral BC (S31) and an auxiliary state
% (S6, g = 1)
S31s = mixed_domain_systems('S31_integral_bc_2d');
S6s = mixed_domain_systems('S6_indomain_1d_2d',struct('g',1));
clear stateNameGenerator
x = pde_var(1,[s1;s2],[0,1;0,1]);
PDEc = initialize([diff(x,t)==diff(x,s1,2)+diff(x,s2,2);
        subs(diff(x,s1),s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0]);
SP = {'Atotal: Neumann-Dirichlet control', PDEc, cos(pi*sx/2)*sin(pi*sy), -(1/4+1)*pi^2;
      'Atotal: S31_integral_bc_2d', initialize(S31s.PDE), cos(pi*sx)*sin(pi*sy), -2*pi^2;
      'Atotal: S6_indomain_1d_2d', initialize(S6s.PDE), [sin(pi*sx), sin(pi*sx)*(1+cos(pi*sy)/2)], -pi^2+1};
for i=1:size(SP,1)
    r = struct('part','piesim','name',SP{i,1},'status','','detail','');
    try
        lam = piesim_spectrum(SP{i,2},mk(16),struct('ic',SP{i,3}));
        r.status = ternary(abs(lam-SP{i,4})<1e-6*abs(SP{i,4}),'pass','FAIL');
        r.detail = sprintf('rightmost Re eig %+.6g, exact %+.6g',lam,SP{i,4});
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
% Systems of the catalogue against the collocation model; initial
% conditions satisfy the BCs, PDE declaration order
opts = mk(16);
L = {'S3_edge1d_heat2d',struct(),{sin(pi*sx), sin(pi*sx)*(1-sy)+sin(pi*sx)*sin(2*pi*sy)/2};
     'S4_wentzell',struct('c',1,'r',0),{sin(pi*sx), sin(pi*sx)*(1-sy)+sin(pi*sx)*sin(2*pi*sy)/2};
     'S5_chain_ode_1d_2d',struct('a',1),{sym(1), (sx-1)+cos(pi*sx/2)/2, (1-sy)*((sx-1)+cos(pi*sx/2)/2)+cos(pi*sx/2)*sin(pi*sy)/2};
     'S8_two_edges',struct('r3',0),{sin(pi*sx), sin(pi*sy), (1-sy)*sin(pi*sx)+sx*sin(pi*sy)+3*sin(pi*sx)*sin(pi*sy)/10};
     'S6_indomain_1d_2d',struct('g',1),{sin(pi*sx), sin(pi*sx)*(1+cos(pi*sy)/2)};
     'S9_ode_2d_flux',struct('a',0.5,'c',1),{sym(1), (1-sy)+cos(pi*sx)*sin(pi*sy)/2};
     'S2_ode_heat2d_edge',struct('a',1),{sym(1), (1-sy)+cos(pi*sx)*sin(pi*sy)/2}};
for i=1:size(L,1)
    r = struct('part','piesim','name',L{i,1},'status','','detail','');
    try
        S = mixed_domain_systems(L{i,1},L{i,2});
        uinput = struct('ic',[L{i,3}{:}]);
        [sol,grd] = PIESIM(initialize(S.PDE),opts,uinput);
        err = err_colloc(sol,grd,L{i,1},L{i,2},L{i,3},tf);
        r.status = ternary(max(err)<tol,'pass','FAIL');
        r.detail = ['errors ',mat2str(err',3)];
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
% An input and an output: S25 (the Wentzell edge with w(t,s1) on x1,
% z = int x1), w = sin(pi s1) sin(2t); reference: the collocation model of
% S4 driven through the x1 rows, the forcing generated by the oscillator
% p' = 2q, q' = -2p (p = sin 2t), so that expm is exact in time
r = struct('part','piesim','name','S25_hinf_wentzell, input and output','status','','detail','');
try
    S = mixed_domain_systems('S25_hinf_wentzell',struct('c',1,'r',0));
    ic = {sin(pi*sx), sin(pi*sx)*(1-sy)+sin(pi*sx)*sin(2*pi*sy)/2};
    uinput = struct('ic',[ic{:}]);  uinput.w = sin(pi*sx)*sin(2*st);
    [sol,grd] = PIESIM(initialize(S.PDE),opts,uinput);
    [err,ez] = err_colloc_forced(sol,grd,ic,tf);
    r.status = ternary(max([err; ez])<tol,'pass','FAIL');
    r.detail = sprintf('states %s, output z %.2e',mat2str(err',3),ez);
catch ME
    r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
end
R(end+1) = report(r);
% Library PDE 43, both signs: mode sin(s1) of x1 obeys
% a' = 0.9 a + sgn a(t-1); the history h(t) = 1 + 0.3 t + g t^2 + k t^3 on
% [-1,0] satisfies h'(0) = 0.9 + sgn h(-1) and h''(0) = 0.27 + sgn h'(-1),
% (g,k) = (-0.11,-0.01) for sgn -1 and (4.47,5.77) for sgn +1, so the
% solution is twice continuously differentiable at t = 0 (a kink would
% travel through the transport state and limit the spectral accuracy)
for sgn=[-1 1]
    r = struct('part','piesim','name',sprintf('S18_lib43_interior_delay sgn %+d',sgn),'status','','detail','');
    try
        if sgn<0, g = -0.11; k = -0.01; else, g = 4.47; k = 5.77; end
        S = mixed_domain_systems('S18_lib43_interior_delay',struct('tau',1,'sgn',sgn));
        uinput = struct('ic',[sin(sx), sin(sx)*(1-3*sy/10+g*sy^2-k*sy^3)]);
        [sol,grd] = PIESIM(initialize(S.PDE),opts,uinput);
        so = dde23(@(t,a,Z) 0.9*a+sgn*Z,1,@(t) 1+0.3*t+g*t^2+k*t^3,[0 tf],ddeset('RelTol',1e-12,'AbsTol',1e-14));
        aT = deval(so,tf);
        err = max(abs(sol.final.primary{2}(:,1)-aT*sin(grd{1}(:))))/abs(aT);
        r.status = ternary(err<tol,'pass','FAIL');
        r.detail = sprintf('error %.2e, a(tf) = %.10g',err,aT);
    catch ME
        r.status = 'FAIL';  r.detail = sprintf('%s (%s:%d)',ME.message,ME.stack(1).name,ME.stack(1).line);
    end
    R(end+1) = report(r);                                                  %#ok<AGROW>
end
end

function lam = piesim_spectrum(PDE,opts,uinput)
% LAM = PIESIM_SPECTRUM(PDE,OPTS,UINPUT) returns the largest real part of
% the eigenvalues of the discrete generator Atotal = Tcheb^{-1} Acheb that
% PIESIM integrates, for a 2-D 'pde_struct' PDE: the steps of PIESIM.m up
% to PIESIM_discretize_all_2D.
[opts,uinput] = PIESIM_options_check({PDE,opts,uinput});
[PDE,uinput,opts,psize] = PIESIM_input_check(PDE,uinput,opts);
PIE = convert_PIETOOLS_PDE(PDE,{'pie','silent'});
if isa(PIE,'pie_struct'), PIE = pie2piesim(PIE); end
PIE = rescalePIE_2D(PIE,[-1,1;-1,1]);
[uinput,psize] = PIESIM_initial_setup(uinput,psize,PIE,opts.type);
Dop = PIESIM_discretize_all_2D(PIE,uinput,psize);
lam = max(real(eig(Dop.Atotal)));
end

function kinds = state_kinds(PDE)
% KINDS = STATE_KINDS(PDE) returns 'ode', 'x', 'y' or 'xy' for each state
% of the 2-D 'pde_struct' PDE, in declaration order.
kinds = cell(numel(PDE.x),1);
for k=1:numel(PDE.x)
    dep = PDE.x_tab(k,3:4);
    if ~any(dep),               kinds{k} = 'ode';
    elseif dep(1) && ~dep(2),   kinds{k} = 'x';
    elseif ~dep(1) && dep(2),   kinds{k} = 'y';
    else,                       kinds{k} = 'xy';
    end
end
end

function err = err_exact(sol,grd,ex,kinds,tf)
% ERR = ERR_EXACT(SOL,GRD,EX,KINDS,TF) returns, per state, the largest
% error of the PIESIM final state on its grid GRD relative to the largest
% value of the exact solution EX (symbolic in st, sx, sy) at time TF.
syms st sx sy real
[X,Y] = ndgrid(grd{1},grd{2});
cnt = struct('ode',0,'x',0,'y',0,'xy',0);
err = zeros(numel(kinds),1);
for k=1:numel(kinds)
    kd = kinds{k};  cnt.(kd) = cnt.(kd)+1;  j = cnt.(kd);
    f = matlabFunction(subs(ex(k),st,tf)+0*sx+0*sy,'Vars',[sx sy]);
    switch kd
        case 'ode', num = sol.final.primary{1}(j);       ref = f(0,0);
        case 'x',   num = sol.final.primary{2}(:,j);     ref = f(grd{1}(:),0*grd{1}(:));
        case 'y',   num = sol.final.primary{3}(:,j);     ref = f(0*grd{2}(:),grd{2}(:));
        case 'xy',  num = sol.final.primary{4}(:,:,j);   ref = f(X,Y);
    end
    err(k) = max(abs(num(:)-ref(:)))/max(1e-14,max(abs(ref(:))));
end
end

function err = err_colloc(sol,grd,name,par,ic,tf)
% ERR = ERR_COLLOC(SOL,GRD,NAME,PAR,IC,TF) returns, per state, the largest
% error of the PIESIM final state relative to the collocation solution at
% time TF from the initial conditions IC (cell of symbolic expressions in
% sx, sy), interpolated to the PIESIM grid.
syms sx sy real
M = mixed_domain_colloc(name,par,28);
x = M.nodes(:);
u0 = zeros(numel(M.idxI)+numel(M.idxB),1);
for k=1:numel(M.state)
    f = matlabFunction(ic{k}+0*sx+0*sy,'Vars',[sx sy]);
    switch M.state(k).kind
        case 'ode', v = f(0,0);
        case 's1',  v = f(x,0*x);
        case 's2',  v = f(0*x,x);
        case 's12', [X1,X2] = ndgrid(x,x);  v = f(X1,X2);
    end
    u0(M.state(k).idx) = v(:);
end
uI = expm(M.L*tf)*u0(M.idxI);
u = zeros(size(u0));  u(M.idxI) = uI;  u(M.idxB) = M.G*uI;
cnt = struct('ode',0,'s1',0,'s2',0,'s12',0);
err = zeros(numel(M.state),1);
for k=1:numel(M.state)
    kd = M.state(k).kind;  cnt.(kd) = cnt.(kd)+1;  j = cnt.(kd);
    v = u(M.state(k).idx);
    switch kd
        case 'ode', num = sol.final.primary{1}(j);     rv = v;
        case 's1',  num = sol.final.primary{2}(:,j);   rv = bary(x,grd{1}(:))*v;
        case 's2',  num = sol.final.primary{3}(:,j);   rv = bary(x,grd{2}(:))*v;
        case 's12'
            num = sol.final.primary{4}(:,:,j);
            rv = bary(x,grd{1}(:))*reshape(v,numel(x),numel(x))*bary(x,grd{2}(:))';
    end
    err(k) = max(abs(num(:)-rv(:)))/max(1e-14,max(abs(rv(:))));
end
end

function [err,ez] = err_colloc_forced(sol,grd,ic,tf)
% [ERR,EZ] = ERR_COLLOC_FORCED(SOL,GRD,IC,TF) compares the PIESIM solution
% of S25 (c = 1, r = 0, w = sin(pi s1) sin(2t), z = int x1) at time TF with
% the collocation model of S4 driven by w through the rows of x1: per state
% the largest error relative to the largest reference value, and the
% relative error of the regulated output z.
syms sx sy real
M = mixed_domain_colloc('S4_wentzell',struct('c',1,'r',0),28);
x = M.nodes(:);  n = numel(x);  nf = numel(M.idxI)+numel(M.idxB);
u0 = zeros(nf,1);
f1 = matlabFunction(ic{1}+0*sx+0*sy,'Vars',[sx sy]);
f2 = matlabFunction(ic{2}+0*sx+0*sy,'Vars',[sx sy]);
u0(M.state(1).idx) = f1(x,0*x);
[X1,X2] = ndgrid(x,x);  v = f2(X1,X2);  u0(M.state(2).idx) = v(:);
% w enters the x1 equation at the interior nodes of x1
bf = zeros(nf,1);  i1 = M.state(1).idx(2:n-1);  bf(i1) = sin(pi*x(2:n-1));
B = bf(M.idxI);  m = numel(M.idxI);
Aug = [M.L, B, zeros(m,1); zeros(1,m), 0, 2; zeros(1,m), -2, 0];
y = expm(Aug*tf)*[u0(M.idxI); 0; 1];
u = zeros(nf,1);  u(M.idxI) = y(1:m);  u(M.idxB) = M.G*y(1:m);
x1 = u(M.state(1).idx);  x2 = reshape(u(M.state(2).idx),n,n);
r1 = bary(x,grd{1}(:))*x1;
r2 = bary(x,grd{1}(:))*x2*bary(x,grd{2}(:))';
n1 = sol.final.primary{2}(:,1);  n2 = sol.final.primary{4}(:,:,1);
err = [max(abs(n1-r1))/max(abs(r1)); max(abs(n2(:)-r2(:)))/max(abs(r2(:)))];
zr = M.w(:)'*x1;
ez = abs(sol.final.regulated{1}(1)-zr)/abs(zr);
end

function P = bary(x,t)
% P = BARY(X,T) returns the barycentric interpolation matrix from the
% Chebyshev-Lobatto nodes X to the points T.
n = numel(x);  w = (-1).^(0:n-1)';  w([1 n]) = w([1 n])/2;
P = zeros(numel(t),n);
for k=1:numel(t)
    d = t(k)-x;  hit = find(abs(d)<1e-14,1);
    if ~isempty(hit), P(k,hit) = 1; continue, end
    q = w./d;  P(k,:) = q'/sum(q);
end
end



%% ------------------------------------------------------------------------
%% Shared helpers
%% ------------------------------------------------------------------------
function r = report(r)
% R = REPORT(R) prints one result line and returns R.
fprintf('%-8s %-34s %-8s %s\n',r.part,r.name,r.status,r.detail);
end

function s = ternary(c,a,b)
% S = TERNARY(C,A,B) returns A if C, else B.
if c, s = a; else, s = b; end
end

function tf = nonempty(P)
% TF = NONEMPTY(P) is true for a 'copvar' with blocks.
tf = ~isempty(P) && isa(P,'copvar') && ~isempty(P.C);
end

function y = apply_cop(P,v,M)
% Y = APPLY_COP(P,V,M) applies the 'copvar' P to the cell V of 'polynomial'
% blocks, each block by its kernels (apply_sopvar): y{i} = sum_j P.C{i,j}
% v{j}. An absent operator gives M scalar zeros.
if ~nonempty(P)
    if nargin<3, M = 0; end
    y = repmat({polynomial(0)},M,1);
    return
end
y = cell(size(P.C,1),1);
for i=1:size(P.C,1)
    acc = polynomial(zeros(P.dim_out(i),1));
    for j=1:size(P.C,2)
        if isempty(P.C{i,j}), continue, end
        acc = acc+apply_sopvar(P.C{i,j},polynomial(v{j}));
    end
    y{i} = acc;
end
end

function c = addc(a,b)
% C = ADDC(A,B) adds two cells of 'polynomial' elementwise.
c = a;
for i=1:numel(a), c{i} = a{i}+b{i}; end
end

function y = Dop(p,pvars,d)
% Y = DOP(P,PVARS,D) differentiates P D(j) times in the variable PVARS{j}.
y = polynomial(p);
for j=1:numel(d)
    for i=1:d(j)
        y = diff(y,polynomial(pvars(j)));
    end
end
end

function y = Fop(p,drop,vreg,dom)
% Y = FOP(P,DROP,VREG,DOM) integrates P over the variables DROP, each on its
% interval (DOM rows ordered as VREG).
y = polynomial(p);
for j=1:numel(drop)
    kv = find(strcmp(vreg,drop{j}));
    y = int(y,polynomial(drop(j)),dom(kv,1),dom(kv,2));
end
end

function e = rel(a,b)
% E = REL(A,B) returns max|coef(A-B)| / max(1,max|coef(B)|).
e = pmax(polynomial(a)-polynomial(b))/max(1,pmax(b));
end

function e = cmpc(a,b)
% E = CMPC(A,B) returns the largest coefficient of A{j}-B{j} relative to
% the largest of B.
sc = 1;  e = 0;
for j=1:numel(a), sc = max(sc,pmax(b{j})); end
for j=1:numel(a), e = max(e,pmax(polynomial(a{j})-polynomial(b{j}))); end
e = e/sc;
end

function m = pmax(p)
% M = PMAX(P) returns the largest absolute coefficient of the 'polynomial' P.
p = polynomial(p);  c = p.coefficient;
if isempty(c), m = 0; else, m = full(max(abs(c(:)))); end
end

function m = maxc(X)
% M = MAXC(X) returns the largest PMAX over the cell X.
m = 0;
for j=1:numel(X), m = max(m,pmax(X{j})); end
end

function p = randpoly(vars,deg,m)
% P = RANDPOLY(VARS,DEG,M) returns an M x 1 'polynomial' in the variables
% VARS (cellstr), of degree at most DEG in each, with normal coefficients.
mons = {polynomial(1)};
for j=1:numel(vars)
    v = polynomial(vars(j));  new = {};
    for i=1:numel(mons)
        for q=0:deg, new{end+1} = mons{i}*v^q; end                         %#ok<AGROW>
    end
    mons = new;
end
p = polynomial(zeros(m,1));
for i=1:numel(mons), p = p+randn(m,1)*mons{i}; end
end

function f = rand_blocks(P,side)
% F = RAND_BLOCKS(P,SIDE) returns random 'polynomial' blocks on the input
% ('in') or output ('out') spaces of the 'copvar' P.
if strcmp(side,'in'), sp = P.space_in; dm = P.dim_in; else, sp = P.space_out; dm = P.dim_out; end
f = cell(size(sp,1),1);
for b=1:size(sp,1), f{b} = randpoly(P.vars(sp(b,:)),2,dm(b)); end
end

function v = ip(u,g,P,side)
% V = IP(U,G,P,SIDE) returns sum_b int u_b' g_b over the variables of
% component b of the input or output space of P.
if strcmp(side,'in'), sp = P.space_in; else, sp = P.space_out; end
v = 0;
for b=1:numel(u)
    q = polynomial(u{b})'*polynomial(g{b});
    for k=find(sp(b,:))
        q = int(q,polynomial(P.vars(k)),P.dom(k,1),P.dom(k,2));
    end
    q = polynomial(q);  c = q.coefficient;
    if ~isempty(c), v = v+full(sum(c(:))); end
end
end
