function R = heatNd_poincare(pie,spec,popts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = HEATND_POINCARE(PIE,SPEC,POPTS) the PIE-form Poincare LPI of a
% HEATND_PIE instance: a diagnostic of the heat stability LPI (HEATND_LPI)
% that shares its T and A but runs SEPARATELY from it and is much smaller.
%
% FORMS. u = T v (v = D^delta u); operators from HEATND_POINCARE_OPS, sums
% over the directions i in SPEC.sel ('full' = all, or one index i):
%   'L2'      Q(lam) = sum_i (D_iT)*(D_iT) - lam T*T >= 0,
%             i.e. ||d_sel u||^2 >= lam ||u||^2 for every u = T v;
%   'energy'  Q(lam) = A_sel*A_sel - lam sum_i (D_iT)*(D_iT) >= 0,
%             i.e. ||sum_i d_i^2 u||^2 >= lam ||d_sel u||^2 (A_sel = sum_i A_i).
% SPEC.target 'grad' (default) writes -T*A_sel as sum_i (D_iT)*(D_iT) (self-
% adjoint by construction; no integration by parts); 'op' uses -(T*A_sel +
% A_sel*T)/2, the same operator by I0 (c) (HEATND_POINCARE_OPS).
%
% EXACT (reference study, DERIVED; both forms): per direction mu_i, the
% lowest eigenvalue of -d^2/ds_i^2 under bc{i} (pi^2/L^2 DD, pi^2/(4L^2)
% DN/ND); full lambda1 = sum_i mu_i. Per direction T*T = T_i^2 (x) T_ihat^2,
% T*A_i = T_i (x) T_ihat^2, A_i*A_i = I (x) T_ihat^2, so the L2 form is
% (-T_i - lam T_i^2) (x) T_ihat^2 and the energy form (I + lam T_i) (x)
% T_ihat^2: >= 0 iff lam <= mu_i. Full: split lam = sum lam_i; the product
% eigenfunction shows lambda1 is sharp. [0,1]^N default: mu_DD = pi^2 =
% 9.869604, mu_DN = pi^2/4 = 2.467401; lambda1 = 9.869604, 12.337006,
% 14.804407 for N = 1, 2, 3.
%
% LINK TO THE HEAT LPI (Cor. 35 at kappa = r + k, HEATND_PIE). P = T gives
% (15a) P*T - ep^2 T*T = (1-ep^2) T*T >= 0 and (15b) P*A + A*P + 2kP*T =
% 2[T*A0 + kappa T*T]: Cor. 35 with P = T certifies EXACTLY kappa <= the
% full L2-form lam. P = -A0 (V = ||grad u||^2) gives (15b) = -2[A0*A0 +
% kappa T*A0]: exactly the full energy form at lam = kappa ((15a) is then
% the L2 form at ep^2 <= lambda1). Both are sharp in exact arithmetic, so a
% heat LPI whose P, R and Q bases contain these certificates certifies at
% least what the Poincare LPI at the matching basis does (README sec. 11).
%
% LPI. Find G = Z*M Z + sum_g Z* g M_g Z, M, M_g >= 0, with G - Q(lam) = 0
% ('lpi_eq_cdopvar' 'symmetric': G and Q(lam) are self-adjoint by
% construction). Z: 'poscopvar' tensor degree SPEC.d (a, c <= d per
% direction, no joint cap), multi-indices
%   L2      all-integral only, {2,3}^N: LOSSLESS IN EXACT ARITHMETIC - the
%           target has no multiplier cell, and a PSD block whose diagonal
%           reaches only zero cells is zero ('eq_opts_sopvar' argument).
%           Not numerically: with the delta indices included, the delta
%           blocks (diagonal ~4e-7, not 0) carry delta x integral cross
%           blocks ~sqrt(tolerance), which inflates certified values ~0.1%
%           and brings the rel_b gate to the edge of certifying above mu
%           (MEASURED, review 09/27/2026: 2-D DDxDN d = 1 faces +1 at 2.002
%           against the pruned 2.0000; 1-D DD d = 3 rel_b 4.9e-7 < 1e-6 at
%           1.00001 mu). Keep this rule for verdict solves; do not widen it;
%   energy  all-integral plus one delta in a direction of sel: A_i*A_i =
%           I (x) T_ihat^2 has the identity in s_i, so all-integral alone
%           would leave its multiplier cell uncovered (0 = b ~= 0).
% Both are the 'eq_opts_sopvar' support rule on the target (R.include_eqopts
% records the comparison). Generators SPEC.gen: 'none'; 'product' ('poscopvar'
% psatz = 1, prod_i (th_i-a_i)(b_i-th_i)); 'faces' (the 2N linear face
% weights (th_i-a_i)/L_i, (b_i-th_i)/L_i by 'poscopvar' psatz = 2i+1 /
% 2i+2 (formerly the copy HEATND_POSW, bitwise the same SDPs), SPEC.faces
% N x 2 logical to restrict; the heat 'linear' Psatz), each its own M_g >= 0
% at the FULL degree d.
%
% LAMBDA. Only b depends on lam (lam multiplies the fixed Q1 in G - Q0 -
% lam Q1). The program is built at three lam (0.5, 0.75, 0.9 of the exact
% value: nonzero, so no row of Q1 is structurally absent), the SDPs are
% checked to have IDENTICAL At, K and m and b affine (third point, <= 1e-12
% relative), and every later lam is b(lam) = b1 + (lam - lam1) db on the
% same At. At is lam-independent, so a HEATND_LINDEP row set (POPTS.lindep)
% is exact for every lam.
%
% MODES (POPTS.mode)
%   'build'     the family and its shape only;
%   'test'      certified verdicts at the absolute lam in POPTS.lams (tests
%               below and above a prediction; no search);
%   'bisect'    largest certified lam (HEATND_BISECT's certified bracket,
%               retry and cap logic, with this family as its solve hook);
%   'objective' one SDP maximizing lam (lam a free variable: At'x - db lam =
%               b(0)), then ONE fixed-lam re-certification at (1 -
%               POPTS.recert) lam-hat. The optimum lies on the boundary of
%               the feasible set, so its own verdict is often uncertain; only
%               the re-certification is a certificate;
%   'sdp'       R.Ds{j}, the SDP at POPTS.lams(j) (HEATND_SOLVE structs
%               sharing one At), and R.D = R.Ds{1}.
% Verdicts are HEATND_SOLVE's (+1 certified feasible, -1 certified
% infeasible - verified Farkas ray, or 0 = b ~= 0 on a row with no variable
% (deficient: the basis misses a target monomial) - 0 uncertain).
%
% INPUT
% - pie:   HEATND_PIE struct (r is irrelevant here: the targets do not
%          involve A's r T part), or an O from HEATND_POINCARE_OPS;
% - spec:  form ('L2'), sel ('full'), d (1), gen ('faces'), faces ([]: all),
%          target ('grad');
% - popts: mode ('build'), lams, solve (HEATND_SOLVE opts, default
%          struct('tight',true): MOSEK tolerances 1e-10 for verdict solves),
%          lindep (false: true solves on independent rows, verdict on all),
%          bisect (HEATND_BISECT bopts overrides; defaults rtol 1e-5, retry
%          {'rows','loose'}), recert (1e-3), cache ('' or a .mat: load the
%          family if present - it must match spec, bc, N and the box dom -
%          else build and save it, -v7.3), keepF (false:
%          true also returns R.F and, for a fresh build, R.aux = program at
%          lams(1), G, Q0, Q1 and RR, for semantic checks of a certificate),
%          verbose (true).
% OUTPUT struct R: spec, N, bc, dom, mu (exact for this sel), lambda1, include,
%   include_eqopts (true if equal to 'eq_opts_sopvar's), sep (its sep),
%   shape (m, nx, Kf, Ks, nnz, zero_rows, uncovered: zero rows of At with b
%   ~= 0 at some lam), affine (At_same, b_err, lams), t (s: ops, target,
%   posvar, eq (per build), sdp, total), mem (MB, process peak private
%   bytes, lifetime), and per mode: test / B (HEATND_BISECT output, lam in
%   the kappa fields: lo, hi, gap = [max(0,mu-hi), mu-lo], trace with MOSEK
%   thread count per attempt), obj (lamhat, v, recert); F, aux if keepF.
%
% Cost (MEASURED shapes, reference study; c = (d+1)(5d+3)): L2 no generator
% m = 2^(N-1) c^N, faces m = 2^(N-1) [c^N + 2N(d+1) c^(N-1)], Gram blocks
% (2(d+1)^2)^N, 1 + 2N with faces. 3-D d = 1 faces: m = 28672, 7 x 512.
% Nothing dense in the decision count: the three builds share one positive
% variable; b updates are O(m). MEASURED 3-D d = 1 (MOSEK 11, 24 threads,
% all rows, tight): all faces ~110 s/iteration (L2) and ~140 s (energy,
% 19.7 GB), i.e. 28-60 min per verdict; faces of one direction only (m
% 20480, 3 x 512) 25-27 s/iteration, 10-11 min per verdict, 9.7 GB; energy
% (m 22528, 3 x 640) 855 s. HEATND_LINDEP (POPTS.lindep) on the 3-D faces At
% exceeded 70 GB: use it in 1-D/2-D only.
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (final reviews): the cache check also compares the box
%   (a family built on another box loaded silently); R.dom stored; the L2
%   'lossless' claim qualified as exact-arithmetic only.
% MMP, 09/27/2026 (library face codes): 'faces' now uses poscopvar psatz =
%   2i+1 / 2i+2, the option copquadvar gained today, not the copy
%   heatNd_posw (retired). MEASURED: the families of the switch-over set
%   (N = 1, 2, 3-D L2 d = 0, 1) bit-identical before/after; caches built
%   before stay valid.
% MMP, 10/01/2026: The program is lpiprogram for every N; it no longer refuses
%   N > 2, so the hand-built copy for N > 2 is commented out. Same program
%   (lpiprogram builds exactly what the copy built).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(spec),   spec = struct();    end
if nargin<3 || isempty(popts),  popts = struct();   end
sd = struct('form','L2','sel','full','d',1,'gen','faces','faces',[],'target','grad');
spec = fillf(spec,sd);
pd = struct('mode','build','lams',[],'solve',struct('tight',true),'lindep',false, ...
            'bisect',struct(),'recert',1e-3,'cache','','keepF',false,'verbose',true);
popts = fillf(popts,pd);
T0 = tic;
if isfield(pie,'DT'),   O = pie;    else,   O = heatNd_poincare_ops(pie);   end
N = O.N;
if ischar(spec.sel) && strcmp(spec.sel,'full'),     I = 1:N;
elseif isnumeric(spec.sel) && isscalar(spec.sel) && any(spec.sel==1:N),  I = spec.sel;
else,   error('heatNd_poincare:sel','SPEC.sel is ''full'' or a direction index 1..%d.',N)
end
mu = sum(O.mu(I));
R = struct('spec',spec,'N',N,'bc',{O.bc},'mu',mu,'lambda1',O.lambda1);

% % % The lam-affine family (built, or loaded from the cache).
if ~isempty(popts.cache) && exist(popts.cache,'file')
    S = load(popts.cache,'F','Rb');
    if ~isequal(S.Rb.spec,spec) || ~isequal(S.Rb.bc,O.bc) || S.Rb.N~=N || ...
       ~isfield(S.Rb,'dom') || ~isequal(S.Rb.dom,O.dom)     % no dom: box unverifiable
        error('heatNd_poincare:cache','Cache %s was built for another spec, PIE or box.',popts.cache)
    end
    F = S.F;    Rb = S.Rb;
    if popts.verbose,   fprintf('heatNd_poincare: family loaded from %s\n',popts.cache);   end
else
    [F,Rb,aux] = build_family(O,spec,I,mu,popts.verbose,popts.keepF);
    if popts.lindep,    F = add_keep(F,popts.verbose);  end
    if ~isempty(popts.cache),   save(popts.cache,'F','Rb','-v7.3');     end
end
if popts.lindep && isempty(F.keep)
    F = add_keep(F,popts.verbose);
    if ~isempty(popts.cache),   save(popts.cache,'F','Rb','-v7.3');     end
end
fn = fieldnames(Rb);
for k = 1:numel(fn),    R.(fn{k}) = Rb.(fn{k});     end

% % % Solves.
so = popts.solve;
if popts.lindep,    so.rows = 'keep';   end
pp = struct('r',0,'exact',struct('lambda1',mu,'kstar',mu));  % HEATND_BISECT: kappa := lam
bo = fillf(popts.bisect,struct('rtol',1e-5,'retry',{{'rows','loose'}},'verbose',popts.verbose));
bo.solve = so;  bo.oracle = @oracle;
switch popts.mode
    case 'build'
    case 'test'
        bo.klist = popts.lams;
        R.test = heatNd_bisect(pp,[],[],[],bo);
    case 'bisect'
        R.B = heatNd_bisect(pp,[],[],[],bo);
    case 'objective'
        R.obj = objective(F,so,popts.recert,popts.verbose);
    case 'sdp'                                  % Ds share At (copy-on-write)
        R.Ds = arrayfun(@(l) sdp_at(F,l),popts.lams,'UniformOutput',false);
        R.D = R.Ds{1};
    otherwise
        error('heatNd_poincare:mode','POPTS.mode is ''build'', ''test'', ''bisect'', ''objective'' or ''sdp''.')
end
R.t.total = toc(T0);    R.mem.peak_private = peakMB();
if popts.keepF,     R.F = F;    if exist('aux','var'),  R.aux = aux;   end,    end

    function v = oracle(lam,s)
        % HEATND_BISECT's solve hook: the certified verdict at lam. Its
        % 'rows' retry asks for the independent rows without supplying them;
        % they are computed once here (exact for every lam: At is fixed).
        if isfield(s,'rows') && strcmp(s.rows,'keep')
            if isempty(F.keep),     F = add_keep(F,popts.verbose);  end
            s.keep = F.keep;
        end
        v = heatNd_solve(sdp_at(F,lam),s);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F,Rb,aux] = build_family(O,spec,I,mu,verbose,keepaux)
% Positive variable once; the equality at three lam; the affine checks.
% AUX (KEEPAUX): the program at lams(1) and the operators, for tests.
aux = [];
N = O.N;    vars = O.vars;  dom = O.dom;
t = struct('ops',0,'target',0,'posvar',0,'eq',[],'sdp',[],'total',0);
% % % Targets Q(lam) = Q0 + lam Q1.
t0 = tic;
GG = [];    As = [];
for i = I
    G = O.DT{i}'*O.DT{i};
    if isempty(GG),     GG = G;     As = O.A{i};    else,   GG = GG + G;    As = As + O.A{i};   end
end
switch spec.target
    case 'grad',    TAs = -1*GG;                            % T*A_sel = -sum (D_iT)*(D_iT)
    case 'op',      TA = O.T'*As;   TAs = 0.5*(TA + TA');   % self-adjoint part
    otherwise,      error('heatNd_poincare:target','SPEC.target is ''grad'' or ''op''.')
end
switch spec.form
    case 'L2',      Q0 = -1*TAs;    Q1 = -1*(O.T'*O.T);
    case 'energy',  Q0 = As'*As;    Q1 = TAs;
    otherwise,      error('heatNd_poincare:form','SPEC.form is ''L2'' or ''energy''.')
end
t.target = toc(t0);
% % % Multi-indices: all-integral, plus one delta in a direction of sel for
% the energy form (header); compared with the eq_opts_sopvar rule.
Ib = dec2bin(0:2^N-1,N)-'0'+2;                  % all-integral, entries {2,3}
Iinc = Ib;
if strcmp(spec.form,'energy')
    for i = I                                   % exactly one delta, in direction i
        Ii = Ib;    Ii(:,i) = 1;    Iinc = [Iinc; unique(Ii,'rows')];   %#ok<AGROW>
    end
end
[eo,~] = eq_opts_sopvar(Q0.C{1,1} + Q1.C{1,1});
inc_eq = isequal(sortrows(Iinc),sortrows(eo.include));
% % % Program and positive variable (one; shared by the three builds).
t0 = tic;
prog = make_prog(vars,dom);
po = struct('include',{{Iinc}});
[prog,G] = poscopvar(prog,1,vars,dom,spec.d,po);
nblk = 1;
switch spec.gen
    case 'none'
    case 'product'
        pp = po;    pp.psatz = 1;
        [prog,G2] = poscopvar(prog,1,vars,dom,spec.d,pp);
        G = G + G2;     nblk = 2;
    case 'faces'
        fc = spec.faces;    if isempty(fc),     fc = true(N,2);     end
        % Face (i,e): poscopvar psatz 2i+1 (e=1) / 2i+2 (e=2), weight       % MMP, 09/27/2026
        % (th_i-a_i)/L_i / (b_i-th_i)/L_i at theta; i indexes the SORTED    % MMP, 09/27/2026
        % registry, which vars = s1..sN is (heatNd_pie: N <= 9), as dom     % MMP, 09/27/2026
        % row order assumes too.                                            % MMP, 09/27/2026
        for i = 1:N
%           si = polynomial(vars(i));   L = dom(i,2)-dom(i,1);              % MMP, 09/27/2026 (was)
%           gi = {(si-dom(i,1))/L, (dom(i,2)-si)/L};                        % MMP, 09/27/2026 (was)
            for e = find(fc(i,:))
%               pg = po;    pg.gfun = gi{e};                                % MMP, 09/27/2026 (was)
%               [prog,Gg] = heatNd_posw(prog,1,vars,dom,spec.d,pg);         % MMP, 09/27/2026 (was)
                pg = po;    pg.psatz = 2*i+e;                               % MMP, 09/27/2026
                [prog,Gg] = poscopvar(prog,1,vars,dom,spec.d,pg);           % MMP, 09/27/2026
                G = G + Gg;     nblk = nblk+1;
            end
        end
    otherwise
        error('heatNd_poincare:gen','SPEC.gen is ''none'', ''product'' or ''faces''.')
end
t.posvar = toc(t0);
% % % Three builds: At must be identical, b affine in lam.
lams = [0.5 0.75 0.9]*mu;
D = cell(1,3);
for j = 1:3
    t0 = tic;
    p = lpi_eq_cdopvar(prog,G - (Q0 + lams(j)*Q1),'symmetric');     t.eq(j) = toc(t0);
    t0 = tic;   D{j} = heatNd_sdp(p);   t.sdp(j) = toc(t0);
    if j==1 && keepaux,    aux = struct('prog',p,'G',G,'Q0',Q0,'Q1',Q1,'RR',D{1}.RR,'lam',lams(1));   end
    clear p
end
same = D{1}.m==D{2}.m && D{1}.m==D{3}.m && isequal(D{1}.K,D{2}.K) && isequal(D{1}.K,D{3}.K) ...
       && isequal(D{1}.At,D{2}.At) && isequal(D{1}.At,D{3}.At);
if ~same,   error('heatNd_poincare:affine','At or K depends on lam: the family is not b-affine.'),   end
bw = cellfun(@(Dj) Dj.b*Dj.bscl,D,'UniformOutput',false);          % raw b
db = (bw{2}-bw{1})/(lams(2)-lams(1));
berr = norm(bw{3} - (bw{1} + (lams(3)-lams(1))*db))/max(norm(bw{3}),eps);
if berr>1e-12,  error('heatNd_poincare:affine','b is not affine in lam (third-point error %.1e).',berr),   end
zr = full(max(abs(D{1}.At),[],1))'<=1e-12;      % sossolve's zero-row test (HEATND_SDP)
unc = zr & (abs(bw{1})>1e-12*max(abs(bw{1})) | abs(db)>1e-12*max(abs(db)));
F = struct('At',D{1}.At,'c',D{1}.c,'K',D{1}.K,'m',D{1}.m,'nx',D{1}.nx,'nnz',D{1}.nnz, ...
           'b1',bw{1},'db',db,'lam1',lams(1),'zr',zr,'keep',[]);
Rb = struct('spec',spec,'N',N,'bc',{O.bc},'dom',dom,'include',Iinc,'include_eqopts',inc_eq,'sep',eo.sep, ...
            'nblocks',nblk, ...
            'shape',struct('m',F.m,'nx',F.nx,'Kf',F.K.f,'Ks',F.K.s,'nnz',F.nnz, ...
                           'zero_rows',nnz(zr),'uncovered',nnz(unc)), ...
            'affine',struct('At_same',same,'b_err',berr,'lams',lams),'t',t, ...
            'mem',struct('peak_private',peakMB()));
if verbose
    fprintf(['heatNd_poincare: N=%d %s sel=%s d=%d gen=%s: m %d, Ks %s, Kf %d, nnz %d, uncovered %d; ' ...
             'At identical at 3 lam, b affine err %.1e; include = eq_opts %d; posvar %.1f s, eq %.1f s\n'], ...
            N,spec.form,num2str(spec.sel),spec.d,spec.gen,F.m,mat2str(F.K.s),F.K.f,F.nnz,nnz(unc), ...
            berr,inc_eq,t.posvar,mean(t.eq));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = sdp_at(F,lam)
% The SDP at lam: b = b1 + (lam - lam1) db, normalized as HEATND_SDP does
% (b/||b||); deficiency recomputed, since it may depend on lam (a T*T
% monomial the basis misses is 0 = lam c).
bw = F.b1 + (lam-F.lam1)*F.db;
s = norm(bw);   if s==0 || ~isfinite(s),    s = 1;  end
D = struct('At',F.At,'b',bw/s,'c',F.c,'K',F.K,'bscl',s,'m',F.m,'nx',F.nx,'nnz',F.nnz, ...
           'deficient',full(any(F.zr & abs(bw)>1e-12*max(abs(bw)))),'lam',lam);
end


function F = add_keep(F,verbose)
% Independent rows of the lam-independent At (HEATND_LINDEP), once.
[F.keep,li] = heatNd_lindep(F.At);
if verbose,     fprintf('  lindep: rank %d of %d rows (%.1f s)\n',li.rank,li.m,li.t);  end
end


function out = objective(F,so,recert,verbose)
% max lam s.t. At'x - db lam = b(0), x in K: lam is a new FIRST free
% variable (SeDuMi order: K.f block first), its column scaled with b so it
% is read off unscaled; then one fixed-lam re-certification.
b0 = F.b1 - F.lam1*F.db;    s = norm(b0);   if s==0,    s = 1;  end
K = F.K;    K.f = K.f + 1;
zo = F.zr & abs(F.db)<=1e-12*max(abs(F.db));    % still no variable, lam included
D = struct('At',[-(F.db(:)/s)'; F.At],'b',b0/s,'c',[-1; zeros(F.nx,1)],'K',K,'bscl',s, ...
           'm',F.m,'nx',F.nx+1,'nnz',nnz(F.At)+nnz(F.db), ...
           'deficient',full(any(zo & abs(b0)>1e-12*max(abs(b0)))));
so.keepx = true;    so.rows = 'all';    so.keep = [];
v = heatNd_solve(D,so);
lamhat = NaN;   if ~isempty(v.x),   lamhat = v.x(1);    end
v.x = [];
out = struct('lamhat',lamhat,'v',v,'recert',[],'lam_recert',NaN);
if isfinite(lamhat) && lamhat>0
    so.keepx = false;
    out.lam_recert = (1-recert)*lamhat;
    out.recert = heatNd_solve(sdp_at(F,out.lam_recert),so);
end
if verbose
    rs = NaN;   if ~isempty(out.recert),    rs = out.recert.st;     end
    fprintf('  objective: lam-hat %.8g (st %+d, %s); re-certification at %.8g: st %+d\n', ...
            lamhat,v.st,v.why,out.lam_recert,rs);
end
end


function prog = make_prog(vars,dom)
% LPI program in the N variables. lpiprogram refuses N > 2 ('more than 2
% spatial variables'); for N > 2 build what it would return, as HEATND_LPI.
% (10/01/2026: it no longer refuses; one call for every N.)                 % MMP, 10/01/2026
% if numel(vars)<=2                                                         % MMP, 10/01/2026 (was)
    prog = lpiprogram(polynomial(vars(:)),[],dom);
% else                                                                      % MMP, 10/01/2026 (was)
%     prog = sosprogram(polynomial([]),dpvar(zeros(0,1)));                  % MMP, 10/01/2026 (was)
%     prog.vartable = [prog.vartable; polynomial(vars(:)); polynomial(strcat(vars(:),'_dum'))]; % MMP, 10/01/2026 (was)
%     prog.dom = dom;                                                       % MMP, 10/01/2026 (was)
% end                                                                       % MMP, 10/01/2026 (was)
end


function s = fillf(s,d)
% Fields of D missing from S.
fn = fieldnames(d);
for i = 1:numel(fn),    if ~isfield(s,fn{i}),   s.(fn{i}) = d.(fn{i});  end,    end
end


function mb = peakMB()
mb = NaN;
try
    p = System.Diagnostics.Process.GetCurrentProcess();     p.Refresh();
    mb = double(p.PeakPagedMemorySize64)/2^20;
catch
end
end
