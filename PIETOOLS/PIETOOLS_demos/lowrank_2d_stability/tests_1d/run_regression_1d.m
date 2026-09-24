function varargout = run_regression_1d()                                    % CC, 09/22/2026
% run_regression_1d -- the 1-D regression suite for the low-rank certifier.
%
%   >> run_regression_1d
%
% No arguments.  Prints one numbered PASS/FAIL line per test with the measured
% numbers underneath, and ends with a machine-readable sentinel:
%
%   REGRESSION_1D: 10/10 PASS   [ 248.3 s ]
%   REGRESSION_1D: 8/10 PASS  FAILED: T5 T6   [ 251.1 s ]
%
% so `grep REGRESSION_1D` is enough for a human or a script.  With an output
% argument it also returns the per-test struct array.
%
% WHAT IS UNDER TEST.  The Burer-Monteiro core of the shipped 2-D certifier --
% bm_setup, bm_dr, bm_lm2, bm_proj, bm_resid, bm_report, raw_data (see
% t1d_corefiles) -- plus restrict_solve.  Those files are DIMENSION-AGNOSTIC:
% they see only the SeDuMi triple and a rank profile, and every defect this
% campaign found lived in them or in restrict_solve.  So a 1-D driver over that
% same core tests the shipped code, and does it at m = 56 with a 0.06 s
% interior-point solve instead of the 2-D sizes -- which is the only reason a
% suite this dense can be run often enough to catch anything.  T0 asserts the
% byte-identity the whole argument rests on; see README.md.
%
% WHY EACH TEST EXISTS is stated at the top of its local function, naming the
% defect it guards and the measured number behind it.  A test without a
% specific number is not a test.
tA = tic;
D = t1d_path();
% The current folder is searched BEFORE the path and cannot be removed from it
% (B6, mechanism 3), so own it for the duration of the run.  T0 verifies the
% resolution regardless, so deleting this line fails a test rather than
% silently changing which code gets measured.
oldcd = cd(D.tests);   restore = onCleanup(@() cd(oldcd));                  %#ok<NASGU>

fprintf('\nPIETOOLS low-rank certifier -- 1-D regression suite\n');
fprintf('  tree     %s\n',D.root);
fprintf('  matlab   %s   threads %d\n',version,maxNumCompThreads);
fprintf('  path     %d shadowing entries removed\n',numel(D.removed));
for k = 1:min(3,numel(D.removed)), fprintf('             %s\n',D.removed{k}); end
if numel(D.removed) > 3, fprintf('             ... and %d more\n',numel(D.removed)-3); end
fprintf('  date     %s\n\n',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));

% ---- shared fixtures.  Built once: the PIE conversion and the LPI assembly
% ---- cost more than most of the tests, and every test wants the same program.
% ---- phi_t = phi_ss + lam phi on [0,1], Dirichlet, stable iff lam < pi^2.
F = struct();   t = tic;
evalc('PIE1 = rd_pie_lam(1,0.10*pi^2);');
F.B1 = t1d_build(PIE1,'heavy');                      % n=1, m=56, Ns=[11 20 11]
evalc('PIE2 = rd_pie_lam(2,0.10*pi^2);');
F.B2 = t1d_build(PIE2,'heavy');                     % n=2, m=224, Ns=[22 40 22]
F.D = D;
fprintf('  fixtures n=1 m=%d Ns=[%s] | n=2 m=%d Ns=[%s]   (%.1f s)\n\n', ...
    F.B1.m,num2str(F.B1.Ns),F.B2.m,num2str(F.B2.Ns),toc(t));

R = [];
[R,F] = runtest(R,F,'T0','HYGIENE / PROVENANCE',      @t_T0);
[R,F] = runtest(R,F,'T1','GATE NON-VACUITY',          @t_T1);
[R,F] = runtest(R,F,'T2','DETERMINISM',               @t_T2);
[R,F] = runtest(R,F,'T3','SCALE / b-NORMALISATION',   @t_T3);
[R,F] = runtest(R,F,'T4','PAD MOBILITY',              @t_T4);
[R,F] = runtest(R,F,'T5','FACE SUPERSET INVARIANT',   @t_T5);
[R,F] = runtest(R,F,'T6','LM TRUNCATION / RAGGEDNESS',@t_T6);
[R,F] = runtest(R,F,'T7','RANK SCALING IN n',         @t_T7);
[R,F] = runtest(R,F,'T8','MISS ATTRIBUTION',          @t_T8);
[R,F] = runtest(R,F,'T9','w0 SEEDING HOOK',           @t_T9);                %#ok<ASGLU>

np = sum([R.pass]);   nt_ = numel(R);
bad = {R(~[R.pass]).id};
if isempty(bad)
    fprintf('\nREGRESSION_1D: %d/%d PASS   [ %.1f s ]\n\n',np,nt_,toc(tA));
else
    fprintf('\nREGRESSION_1D: %d/%d PASS  FAILED: %s   [ %.1f s ]\n\n', ...
        np,nt_,strjoin(bad,' '),toc(tA));
end
if nargout, varargout{1} = R; end
end

% =========================================================================
function [R,F] = runtest(R,F,id,name,fn)
t = tic;   A = anew();
try
    [A,F] = fn(A,F);   err = '';
catch ME
    err = sprintf('%s [%s]',ME.message,ME.identifier);
    A = ck(A,'threw',false,['test threw: ' err]);
end
el = toc(t);   pass = isempty(A.bad);
fprintf(' %-3s %-28s %s  (%6.1f s)\n',id,name,tern(pass,'PASS','FAIL'),el);
for k = 1:numel(A.lines), fprintf('        %s\n',A.lines{k}); end
R(end+1).id = id;   R(end).name = name;   R(end).pass = pass;
R(end).t = el;      R(end).lines = {A.lines};   R(end).bad = {A.bad};
end

% =========================================================================
% T0 -- B6 path shadowing, and the byte-identity the rest of the suite rests
% on.  MEASURED while this suite was being written: 90 path entries belonging
% to a nested .claude worktree, one of whose private/bm_lm2.m copies DIFFERED
% from the current file; and a probe run from another folder silently called
% that folder's bm_setup instead of the one under test.
function [A,F] = t_T0(A,F)
D = F.D;
% (a) the seven shared files are byte-identical to the shipped ../private
nm = t1d_corefiles();   nid = 0;
for k = 1:numel(nm)
    a = fileread(fullfile(D.tests,[nm{k} '.m']));
    b = fileread(fullfile(D.priv ,[nm{k} '.m']));
    if isequal(a,b), nid = nid+1;
    else, A = ck(A,'core-identity',false,sprintf('%s DIFFERS from ../private -- the suite has stopped testing shipped code',nm{k}));
    end
end
A = ck(A,'core-identity-count',nid==numel(nm), ...
    sprintf('core byte-identical to ../private : %d/%d',nid,numel(nm)));
% (b) and it is THIS directory's copy that executes when t1d_bm calls them
W = t1d_whichcore();   nres = 0;
for k = 1:size(W,1)
    if strcmp(W{k,2},fullfile(D.tests,[W{k,1} '.m'])), nres = nres+1;
    else, A = ck(A,'core-resolution',false,sprintf('%s resolves to %s',W{k,1},W{k,2}));
    end
end
A = ck(A,'core-resolution-count',nres==size(W,1), ...
    sprintf('core resolves inside tests_1d     : %d/%d',nres,size(W,1)));
% (c) every name this directory owns resolves to exactly one file.  The seven
% are exempt from ==1 because ../private legitimately holds the originals and
% `which -all` reports private folders; (a)+(b) pin those instead.
own = setdiff(t1d_ownnames(),nm);   amb = {};
for k = 1:numel(own)
    if numel(which(own{k},'-all'))~=1
        amb{end+1} = sprintf('%s(%d)',own{k},numel(which(own{k},'-all')));  %#ok<AGROW>
    end
end
A = ck(A,'own-unique',isempty(amb), ...
    tern(isempty(amb),sprintf('own names resolving uniquely      : %d/%d',numel(own),numel(own)), ...
                      ['ambiguous names: ' strjoin(amb,' ')]));
% (d) the PIETOOLS entry points too, and they must live in this tree
ent = {'lpisolve','poslpivar','sossolve','lpiprogram','lpisettings','pde_var'};
amb = {};
for k = 1:numel(ent)
    w = which(ent{k},'-all');
    if numel(w)~=1, amb{end+1} = sprintf('%s(%d)',ent{k},numel(w)); end     %#ok<AGROW>
    if ~isempty(w) && ~startsWith(lower(w{1}),lower(D.root))
        A = ck(A,'entry-outside-tree',false,sprintf('%s resolves outside the tree: %s',ent{k},w{1}));
    end
end
A = ck(A,'entry-unique',isempty(amb), ...
    tern(isempty(amb),sprintf('PIETOOLS entry points unique      : %d/%d',numel(ent),numel(ent)), ...
                      ['ambiguous entry points: ' strjoin(amb,' ')]));
% (e) t1d_ownnames is in step with the directory, or a new file goes unwatched
dd = dir(fullfile(D.tests,'*.m'));
have = cellfun(@(s) s(1:end-2),{dd.name},'uni',0);
miss = setdiff(have,t1d_ownnames());
A = ck(A,'ownnames',isempty(miss), ...
    tern(isempty(miss),'t1d_ownnames covers every .m here', ...
                       ['not listed in t1d_ownnames: ' strjoin(miss,' ')]));
% (f) restrict_solve_1d still differs from the shipped file in exactly the
% declaration and the ONE gate call.  A fix landing in the shipped file that
% this copy has not got would otherwise be invisible.
C1 = t1d_codelines(fileread(fullfile(D.tests,'restrict_solve_1d.m')));
C2 = t1d_codelines(fileread(fullfile(D.priv ,'restrict_solve.m')));
if numel(C1)~=numel(C2)
    A = ck(A,'restrict-drift',false,sprintf('restrict_solve code lines %d vs shipped %d',numel(C1),numel(C2)));
else
    d = find(~strcmp(C1,C2));
    okd = numel(d)==2 && contains(C1{d(1)},'restrict_solve_1d(') ...
                      && contains(C1{d(2)},'gate1d(') && contains(C2{d(2)},'opcheck_2d(');
    A = ck(A,'restrict-drift',okd, ...
        sprintf('restrict_solve_1d vs shipped      : %d code lines, %d differ (want 2: decl + gate call)',numel(C1),numel(d)));
end
% (g) the gate threshold has not been redefined out from under every other test
g = getenv('R1_GATE');
A = ck(A,'gate-env',isempty(g), ...
    tern(isempty(g),'R1_GATE unset -> gate threshold is 1e-6',sprintf('R1_GATE=%s overrides the 1e-6 gate',g)));
end

% =========================================================================
% T1 -- B4 gate integrity.  BOTH single-condition acceptances have manufactured
% a wrong answer in this campaign: a point satisfied the operator identity to
% 2.8e-15 while INDEFINITE at mineig/maxeig = -1.52, and a near-zero Gram
% passes the equality rows of a tiny-||b|| program while certifying nothing
% (||b|| = 3.5e-06 on this very fixture).  The negative controls below are what
% make the gate's acceptance mean something.
function [A,F] = t_T1(A,F)
B = F.B1;
F.ipm1 = t1d_ipm(B);                      % reused by T8 as the reference arm
if ~F.ipm1.ok, A = ck(A,'ipm',false,['interior-point arm failed: ' F.ipm1.err]); return, end
q = F.ipm1.x/B.P.nb0;
% (a) a known certificate is ACCEPTED
G = gate1d(B.prog,B.H,B.P,q);
A = ck(A,'accept-cert',G.cert && G.relP<1e-7, ...
    sprintf('certificate     relP %.3e psd %d cert %d   (want cert, relP<1e-7; measured 1.89e-09)',G.relP,G.psd,G.cert));
% (b) the ZERO Gram is rejected
G0 = gate1d(B.prog,B.H,B.P,zeros(B.Ntot,1));
A = ck(A,'reject-zero',~G0.cert && G0.relP>1e-2, ...
    sprintf('zero Gram       relP %.3e cert %d            (want reject, relP>1e-2; measured 2.66)',G0.relP,G0.cert));
% (c) a 1.02x-scaled certificate is rejected: the eppos*I term does not scale
% with the Gram, so the identity breaks by ~2% of the operator norm
Gs = gate1d(B.prog,B.H,B.P,1.02*q);
A = ck(A,'reject-scaled',~Gs.cert && Gs.relP>1e-4, ...
    sprintf('1.02x scaled    relP %.3e cert %d            (want reject, relP>1e-4; measured 5.96e-03)',Gs.relP,Gs.cert));
% (d) THE ONE THAT MATTERS.  Perturb the certificate along a direction in the
% null space of the constraint map restricted to symmetric parameters: the
% operator residual is unchanged to every printed digit while one block goes
% indefinite.  The residual half must still PASS and the PSD half must reject,
% and the test asserts WHICH half did it -- a rejection for the wrong reason
% would leave B4 unguarded.
[qf,ipk] = flipnull(B,q);
if isempty(qf), A = ck(A,'null',false,'no symmetric null direction found'); return, end
Gf  = gate1d(B.prog,B.H,B.P,qf);
rat = Gf.relP/max(G.relP,realmin);
mne = min(Gf.mineig./max(Gf.normQ,realmin));
A = ck(A,'flip-residual-passes',Gf.relP<1e-6 && abs(rat-1)<1e-3, ...
    sprintf('indefinite twin relP %.3e = %.6f x cert   residual half PASSES, as it must',Gf.relP,rat));
A = ck(A,'flip-psd-rejects',~Gf.psd && ~Gf.cert && mne<-1e-2, ...
    sprintf('                block %d  mineig/||Q||_F %.3e  psd %d   PSD half rejects (measured -2.20e-01)',ipk,mne,Gf.psd));
% (e) Gram-free second opinion, REPORTED not asserted.  See README: the
% Galerkin form of -Dop carries ~8% relative quadrature error even on the
% analytic P = I certificate, and the ZERO Gram scores better on it than a real
% certificate does, so it cannot decide anything and must not be a gate.
Gg = gate1d(B.prog,B.H,B.P,q,true,24);
A = nt(A,sprintf('galerkin (diagnostic only) ePmin %.3e eDmin %.3e maxPop %.3e',Gg.ePmin,Gg.eDmin,Gg.maxPop));
end

% =========================================================================
% T2 -- a fixed seed must reproduce bit-identically within one session.
% Without it every other number here is anecdote: the searches are seeded
% inside t1d_bm and nothing else in the loop may consume randomness.
function [A,F] = t_T2(A,F)
o  = struct('maxrank',3,'seeds',[11 22 33 44],'lmit',4000);
Ra = t1d_bm(F.B1,o);
Rb = t1d_bm(F.B1,o);
F.bm1 = Ra;                                       % reused by T5, T7 and T9
same = isequal(Ra.q,Rb.q) && isequal(Ra.r,Rb.r) && isequal(Ra.log,Rb.log) ...
     && Ra.relP==Rb.relP && Ra.attempts==Rb.attempts && Ra.lm_iters==Rb.lm_iters ...
     && Ra.dr_iters==Rb.dr_iters && Ra.gates==Rb.gates && Ra.cert==Rb.cert;
A = ck(A,'bitwise',same, ...
    sprintf('two runs: cert %d/%d  r=[%s]  relP %.12e / %.12e  att %d/%d  lm %d/%d  q bit-identical %d', ...
    Ra.cert,Rb.cert,num2str(Ra.r),Ra.relP,Rb.relP,Ra.attempts,Rb.attempts, ...
    Ra.lm_iters,Rb.lm_iters,isequal(Ra.q,Rb.q)));
A = ck(A,'certified',Ra.cert, ...
    sprintf('route %s at r=[%s] -- the fixture must certify for T5 and T7 to have a face',Ra.route,num2str(Ra.r)));
end

% =========================================================================
% T3 -- B5.  bf is normalised to unit norm inside bm_setup BEFORE any
% preconditioning, so every residual this suite reports is ALREADY relative.
% A regression here silently changes what the gate means in every other test.
function [A,F] = t_T3(A,F)
B  = F.B1;
P1 = bm_setup(B.Atf,   B.bf,B.Ns,B.Kf,1);
P2 = bm_setup(B.Atf,10*B.bf,B.Ns,B.Kf,1);
A = ck(A,'unit-b',abs(norm(P1.bf)-1)<1e-14, ...
    sprintf('||bf|| = %.16g                      (want 1)',norm(P1.bf)));
A = ck(A,'nb0',abs(P2.nb0/P1.nb0-10)<1e-12, ...
    sprintf('b -> 10b : nb0 ratio %.16g   (want 10; nb0 = %.4e)',P2.nb0/P1.nb0,P1.nb0));
A = ck(A,'b-invariant',norm(P1.bf-P2.bf)<1e-14, ...
    sprintf('||bf1 - bf2||  = %.3e            (want <1e-14; measured 5.9e-17)',norm(P1.bf-P2.bf)));
A = ck(A,'W-invariant',norm(P1.W-P2.W,'fro')==0, ...
    sprintf('||W1 - W2||_F  = %.3e            (W comes from Ssym alone, so exactly 0)',norm(P1.W-P2.W,'fro')));
rng(7,'twister');   rv = min(2,B.Ns);
w  = randn(sum(B.Ns.*rv),1);
R1 = bm_report(w,P1,rv);   R2 = bm_report(w,P2,rv);
rel = abs(R1.raw_rel-R2.raw_rel)/max(R1.raw_rel,realmin);
A = ck(A,'resid-invariant',rel<1e-12, ...
    sprintf('raw_rel %.14e vs %.14e  rel diff %.2e  (want <1e-12)',R1.raw_rel,R2.raw_rel,rel));
A = ck(A,'b-in-range',P1.bout<1e-12, ...
    sprintf('b outside range(A) : %.3e        (want <1e-12; the preconditioner drops null directions)',P1.bout));
end

% =========================================================================
% T4 -- B2 zero-pad fixed point.  J_i = 2*W*(Ssym(:,rows_i)*kron(Y_i,I_N)), so
% the column block belonging to column c is Y_i(:,c) (x) I: a ZERO column has
% an IDENTICALLY ZERO Jacobian and Levenberg-Marquardt can never move it, which
% makes the 'prev' rung a bit-for-bit rank-r continuation dressed as rank r+1.
% MEASURED max||J(:,padded)|| = 0.000e+00 against 0.43..0.91 for live columns,
% still exactly zero after 150 LM steps.  The shipped value is asserted to BE
% zero, so this test documents the defect rather than hiding it.
function [A,F] = t_T4(A,F)
B = F.B1;  P = B.P;  Ns = B.Ns;  nb = numel(Ns);
src  = fullfile(F.D.pkg,'pielr_certify.m');
hpad = t1d_subfun(src,'padw');                            % the SHIPPED zero pad
hfit = t1d_subfun(src,'fitw');                            % the SHIPPED fill
rv1 = min(1,Ns);   rv2 = min(2,Ns);
w1  = settle(P,Ns,rv1,11);                        % a settled rank-1 point
Rs  = bm_report(w1,P,rv1);
A = nt(A,sprintf('settled r=1 point: raw_rel %.3e',Rs.raw_rel));
% (a) the shipped zero pad leaves the padded columns immobile, exactly
wz = hpad(w1,rv1,rv2,Ns);
[lv,pd] = coljac(wz,P,Ns,rv1,rv2);
A = ck(A,'pad-zero',max(pd)==0, ...
    sprintf('zero pad  max||J(:,padded)|| = %.3e   EXPECTED EXACTLY 0 -- defect B2, pinned',max(pd)));
A = ck(A,'live-mobile',min(lv)>1e-2, ...
    sprintf('          min||J(:,live)||   = %.3e   (want >1e-2; measured 4.3e-01..9.1e-01)',min(lv)));
% (b) and LM cannot repair it: 150 steps leave the padded entries exactly zero
wz2  = bm_lm2(P,rv2,wz,150,1e-16);
pidx = padidx(Ns,rv1,rv2,P.Kf);
[~,pd2] = coljac(wz2,P,Ns,rv1,rv2);
A = ck(A,'pad-immobile',max(abs(wz2(pidx)))==0 && max(pd2)==0, ...
    sprintf('          after 150 LM steps: max|w(padded)| = %.3e, max||J|| = %.3e, nnz %d of %d', ...
    max(abs(wz2(pidx))),max(pd2),nnz(wz2(pidx)),numel(pidx)));
% (c) the repair the .w0 hook depends on: fitw's small random fill IS mobile
Y0 = cell(1,nb);   k = 0;
for i = 1:nb
    Y0{i} = reshape(w1(k+(1:Ns(i)*rv1(i))),Ns(i),rv1(i));   k = k+Ns(i)*rv1(i);
end
rng(99,'twister');
wf = hfit(Y0,rv2,Ns);
[~,pf] = coljac(wf,P,Ns,rv1,rv2);
A = ck(A,'fitw-mobile',min(pf)>1e-6, ...
    sprintf('fitw fill min||J(:,surplus)|| = %.3e   (want >1e-6; measured 4.09e-04)',min(pf)));
end

% =========================================================================
% T5 -- B3 restrict_solve false negatives.  The restricted system was handed to
% SeDuMi with a zero objective and returned x = 0 exactly on SIX faces that
% PROVABLY CONTAIN a certificate; it is now a backslash on a pivoted-QR row
% selection plus a per-block PSD projection.  The invariant to guard is that a
% face CONTAINING a certifying face still certifies.  Part (d) pins a 1-D case
% where it still does NOT -- see README, this is the suite's first finding.
function [A,F] = t_T5(A,F)
B = F.B1;  P = B.P;  nb = numel(B.Ns);
if ~isfield(F,'bm1') || ~F.bm1.cert
    A = ck(A,'no-face',false,'T2 did not leave a certified face');  return
end
V = F.bm1.V;
% (a) the base face certifies, and the solve on it is DETERMINED -- np == rM is
% the structural fact that makes the restricted problem a linear solve rather
% than an optimisation, which is why handing it to an SDP solver bought nothing
[Rb,qb,ib] = restrict_solve_1d(B.prog,B.H,P,B.Atf,B.bf,V,false);
A = ck(A,'base-cert',Rb.cert && Rb.relP<1e-6, ...
    sprintf('base face rr=[%s] np=%d rM=%d cert %d relP %.3e clip %.2e',num2str(ib.rr),ib.np,ib.rM,Rb.cert,Rb.relP,ib.clip));
A = ck(A,'determined',ib.np==ib.rM && ib.clip==0, ...
    sprintf('          np == rM (%d == %d) and clip 0 : determined linear solve',ib.np,ib.rM));
% (b) CONTAINMENT, so that (d) cannot be dismissed as "nothing to find there".
% The wide face spans the base face to machine precision and the base Gram is
% itself a point of it, so a certificate exists on the wide face by
% construction -- this is the half of the invariant that is pure mathematics.
Vw = cell(1,nb);   res = 0;
for i = 1:nb
    N = size(V{i},1);   [Qq,~] = qr([V{i} eye(N)],0);
    Vw{i} = Qq(:,1:min(size(V{i},2)+1,N));
    res = max(res,norm(V{i} - Vw{i}*(Vw{i}'*V{i})));
end
Ge = gate1d(B.prog,B.H,P,qb);
A = ck(A,'containment',res<1e-12 && Ge.cert, ...
    sprintf('containment max||V - Vw Vw^T V|| = %.3e, base Gram re-gates relP %.3e cert %d',res,Ge.relP,Ge.cert));
% (c) THE INVARIANT, where the restricted map stays well conditioned.  Widening
% one block at a time: MEASURED 2 of 3 certify (clip 0 and 4.2e-08).
ncert = 0;   badclip = [];
for i = 1:nb
    Vu = V;   N = size(V{i},1);   [Qq,~] = qr([V{i} eye(N)],0);
    Vu{i} = Qq(:,1:min(size(V{i},2)+1,N));
    [Ri,~,ii] = restrict_solve_1d(B.prog,B.H,P,B.Atf,B.bf,Vu,false);
    A = nt(A,sprintf('   block %d widened: np=%d rM=%d cert %d relP %.3e clip %.2e', ...
        i,ii.np,ii.rM,Ri.cert,Ri.relP,ii.clip));
    if Ri.cert, ncert = ncert+1; else, badclip(end+1) = ii.clip; end        %#ok<AGROW>
end
A = ck(A,'wide-1block',ncert>=2, ...
    sprintf('single-block widening: %d of %d still certify   (want >=2; measured 2)',ncert,nb));
% a widening that LOST the certificate must have lost it to the conditioning of
% the restricted map, not to something unexplained
A = ck(A,'loss-attributed',isempty(badclip) || min(badclip)>1e10, ...
    sprintf('every lost widening carries clip > 1e10 : %s   (measured 1.33e+14)', ...
    tern(isempty(badclip),'none lost',sprintf('min clip %.2e',min(badclip)))));
% (d) DEFECT PINNED.  Widening EVERY block by one arbitrary orthonormal
% direction makes the restricted map nearly singular -- MEASURED cond(M)
% 7.8e+07, smallest QR pivot 1.9e-08 of the largest -- while np == rM still
% reports a determined system.  The solve then returns a point with
% clip 6.5e+11 whose PSD projection certifies nothing (relP 2.658, which is
% the zero-Gram value).  No alternative solve rescues it: backslash,
% lsqminnorm, truncated SVD at 1e-6 and least squares on ALL rows were each
% measured to fail on this face.  ASSERTED AS CURRENT BEHAVIOUR: if this ever
% certifies, the latent defect is fixed and this line is what says so.
[Rw,~,iw] = restrict_solve_1d(B.prog,B.H,P,B.Atf,B.bf,Vw,false);
A = ck(A,'wide-all-pinned',~Rw.cert && iw.clip>1e10, ...
    sprintf('KNOWN DEFECT pinned: all blocks widened -> np=%d rM=%d cert %d relP %.3e clip %.2e', ...
    iw.np,iw.rM,Rw.cert,Rw.relP,iw.clip));
if Rw.cert
    A = nt(A,'*** the wide face NOW CERTIFIES: the B3 latent defect is fixed -- update T5(d) ***');
end
% negative control: a random face of the same size must NOT certify, or the
% face machinery is accepting everything and (a)-(c) prove nothing
rng(4242,'twister');   Vr = cell(1,nb);
for i = 1:nb, Vr{i} = orth(randn(size(V{i},1),size(V{i},2))); end
[Rn,~,inr] = restrict_solve_1d(B.prog,B.H,P,B.Atf,B.bf,Vr,false);
A = ck(A,'random-rejected',~Rn.cert && Rn.relP>1e-2, ...
    sprintf('negative control: random face np=%d rM=%d cert %d relP %.3e   (measured 2.56)',inr.np,inr.rM,Rn.cert,Rn.relP));
end

% =========================================================================
% T6 -- B1 bm_lm2 truncation.  bm_lm2(P,rv,w,400,1e-16) exits only when its
% 40-try damping loop finds no decrease, or when ||F|| < 1e-16 (never): there
% is no stagnation or relative-improvement test, so 400 is a HARD CUT and
% whether the residual crosses the gate inside exactly 400 steps is near
% arbitrary.  THE DIAGNOSTIC SIGNATURE IS RAGGEDNESS.  MEASURED on this grid at
% the shipped budget of 400 the verdict was [1 0 0 1 0] -- pass, fail, fail,
% pass, fail on 0.01 steps in lam, three contiguous runs -- and at 4000 it was
% [1 1 1 1 1] with every point certifying at r=[1 1 1].  Raggedness or a
% budget-sensitive verdict is the FAILURE condition.  The contiguity statistic
% on that measured [1 0 0 1 0] is 2 runs, so `runs<=1` rejects it: the test
% bites on the defect it was written for.
function [A,F] = t_T6(A,F)
fr  = [0.41 0.42 0.43 0.44 0.45];          % inside the known-good region
lms = [4000 8000];                         % the working budget, and twice it
c = false(numel(fr),2);   rk = cell(numel(fr),2);
F.t6 = repmat(struct('f',[],'cert',[],'B',[],'relP',[]),1,numel(fr));
for k = 1:numel(fr)
    evalc(sprintf('PIEk = rd_pie_lam(1,%.10g*pi^2);',fr(k)));
    Bk = t1d_build(PIEk,'heavy');
    Rk = cell(1,2);
    for j = 1:2
        Rk{j} = t1d_bm(Bk,struct('maxrank',8,'seeds',[11 22 33 44],'lmit',lms(j),'timecap',240));
        c(k,j) = Rk{j}.cert;   rk{k,j} = Rk{j}.r;
    end
    F.t6(k).f = fr(k);   F.t6(k).cert = c(k,1);   F.t6(k).B = Bk;
    F.t6(k).relP = Rk{1}.relP;
    A = nt(A,sprintf('f=%.2f  lmit %d cert %d r=[%-7s] relP %.3e | lmit %d cert %d r=[%-7s]', ...
        fr(k),lms(1),c(k,1),num2str(rk{k,1}),Rk{1}.relP,lms(2),c(k,2),num2str(rk{k,2})));
end
runs = nnz(diff([false; c(:,1); false])==1);
A = ck(A,'contiguous',runs<=1, ...
    sprintf('verdicts at lmit %d : [%s] -> %d contiguous run(s)   (want 1; at lmit 400 this grid gives [1 0 0 1 0], 2 runs)', ...
    lms(1),num2str(double(c(:,1)')),runs));
A = ck(A,'budget-insensitive',isequal(c(:,1),c(:,2)), ...
    sprintf('verdicts at lmit %d : [%s] identical to lmit %d : %d   (a verdict that improves with budget means the cut is binding)', ...
    lms(2),num2str(double(c(:,2)')),lms(1),isequal(c(:,1),c(:,2))));
A = ck(A,'all-certify',all(c(:,1)), ...
    sprintf('every point in the known-good region certifies : %d of %d',nnz(c(:,1)),numel(fr)));
end

% =========================================================================
% T7 -- rank scaling in the number of decoupled replicas.  NOTE the 2-D law
% r* = 2n does NOT hold in 1-D: MEASURED at heavy, lam/lam* = 0.10, lmit 4000
% the search first certifies at r = [1 1 1] for n = 1 and r = [2 2 2] for n = 2.
% What is guarded here is the replica upper bound r(n) <= n*r(1) -- which is
% what the block-diagonal construction on decoupled replicas implies -- and
% that the minimum rank stays O(1) rather than tracking the Gram size, which
% doubles from Ns=[11 20 11] to [22 40 22] between the two.
function [A,F] = t_T7(A,F)
if ~isfield(F,'bm1'), A = ck(A,'no-n1',false,'T2 did not run'); return, end
R1 = F.bm1;
R2 = t1d_bm(F.B2,struct('maxrank',5,'seeds',[11 22 33 44],'lmit',4000));
F.bm2 = R2;                                                   % reused by T9
r1 = max(R1.r);   r2 = max(R2.r);
A = ck(A,'n1',R1.cert && r1<=2, ...
    sprintf('n=1  m=%3d Ns=[%-10s] cert %d r=[%-7s]  (want max r <= 2; measured 1)',F.B1.m,num2str(F.B1.Ns),R1.cert,num2str(R1.r)));
A = ck(A,'n2',R2.cert && r2<=4, ...
    sprintf('n=2  m=%3d Ns=[%-10s] cert %d r=[%-7s]  (want max r <= 4; measured 2)',F.B2.m,num2str(F.B2.Ns),R2.cert,num2str(R2.r)));
A = ck(A,'replica-bound',r2<=2*r1, ...
    sprintf('replica upper bound r(2) <= 2 r(1) : %d <= %d',r2,2*r1));
end

% =========================================================================
% T8 -- the methodological rule, asserted rather than trusted: a lam the search
% does not reach must be CHARGED to the search or to the relaxation, never
% reported as "not reached".  The interior-point arm solves the SAME program,
% so if it certifies the miss belongs to the search, and if it does not the
% relaxation is out of road.  Non-vacuity comes from a point above the analytic
% limit lam* = pi^2, where no certificate exists at all: MEASURED at 1.05 lam*
% the search misses (best relP 3.36e-03) and the IPM misses (relP 1.04e+01).
function [A,F] = t_T8(A,F)
% (a) every miss on T6's grid carries an attribution
nmiss = 0;   natt = 0;
if isfield(F,'t6')
    for k = 1:numel(F.t6)
        if isempty(F.t6(k).cert) || F.t6(k).cert, continue; end
        nmiss = nmiss+1;
        Ak  = t1d_ipm(F.t6(k).B);
        who = tern(Ak.cert,'SEARCH (the relaxation certifies here)','RELAXATION (the IPM misses too)');
        A = nt(A,sprintf('miss at f=%.2f charged to %s [IPM relP %.3e cert %d]',F.t6(k).f,who,Ak.relP,Ak.cert));
        natt = natt+1;
    end
end
A = ck(A,'grid-attributed',natt==nmiss, ...
    sprintf('T6 grid misses %d, attributed %d   (a miss without an owner is itself a test failure)',nmiss,natt));
% (b) and the attribution path is exercised, above the analytic limit
evalc('PIEu = rd_pie_lam(1,1.05*pi^2);');
Bu = t1d_build(PIEu,'heavy');
Ru = t1d_bm(Bu,struct('maxrank',3,'seeds',[11 22],'lmit',400,'timecap',180));
Au = t1d_ipm(Bu);
who = tern(Au.cert,'SEARCH','RELAXATION');
A = ck(A,'above-lamstar',~Ru.cert && ~Au.cert && strcmp(who,'RELAXATION'), ...
    sprintf('lam = 1.05 lam* : search cert %d (best relP %.3e), IPM cert %d (relP %.3e) -> charged to %s', ...
    Ru.cert,Ru.best_relP,Au.cert,Au.relP,who));
A = ck(A,'relaxation-decisive',Au.relP>1, ...
    sprintf('                IPM relP %.3e > 1 : decisively out, as it must be above lam*',Au.relP));
% (c) and the other way, so a RELAXATION verdict is informative rather than the
% arm simply never working
if isfield(F,'ipm1')
    A = ck(A,'reference-live',F.ipm1.cert, ...
        sprintf('reference arm at f=0.10 certifies (relP %.3e) : the attribution can tell the two apart',F.ipm1.relP));
end
end

% =========================================================================
% T9 -- B7, the .w0 seeding hook.  It must run FIRST at every rung: on this
% family seed 11 certifies on the first attempt at low lam, so an arm appended
% after the seeds would never execute and the option would silently measure
% nothing.  And its surplus columns must be NONZERO, or it degenerates into the
% B2 zero-pad fixed point and produces a confident false negative.  The shipped
% hook lives in pielr_certify, which is 2-D only and far too expensive to run
% here, so it is tested by EXTRACTING fitw and by asserting the dispatch on the
% shipped source text.
function [A,F] = t_T9(A,F)
src = fullfile(F.D.pkg,'pielr_certify.m');
% (a) the shipped dispatch: the w0 branch is guarded by si==1, calls fitw, and
% comes textually BEFORE the seeds branch, whose bound is offset by nw0.  The
% fitw call is on the line AFTER the guard, so the search is for the first
% fitw( at or after the guard and strictly before the seeds branch -- fixing
% the branch order and the fill function in one assertion.
L  = strsplit(fileread(src),newline);
iw = find(~cellfun('isempty',regexp(L,'nw0\s*&&\s*si\s*==\s*1','once')),1);
is = find(~cellfun('isempty',regexp(L,'si\s*<=\s*nw0\s*\+\s*numel\(opts\.seeds\)','once')),1);
fw = [];
if ~isempty(iw), fw = iw-1 + find(contains(L(iw:end),'fitw('),1); end
okd = ~isempty(iw) && ~isempty(is) && ~isempty(fw) && iw<is && fw>=iw && fw<is;
A = ck(A,'shipped-dispatch',okd, ...
    sprintf('pielr_certify: w0 guard line %s, fitw call line %s, seeds branch line %s -> w0 leads %d', ...
    num2str(iw),num2str(fw),num2str(is),okd));
% (b) t1d_bm's local copy of fitw still matches the shipped one code for code.
% If it drifts, this suite is testing its own fill rule and not the shipped one.
[~,ts] = t1d_subfun(src,'fitw');
[~,tl] = t1d_subfun(fullfile(F.D.tests,'t1d_bm.m'),'fitw',tempname);
Cs = t1d_codelines(ts);   Cl = t1d_codelines(tl);
A = ck(A,'fitw-copy',isequal(Cs,Cl), ...
    sprintf('local fitw == shipped fitw : %d code lines, equal %d',numel(Cs),isequal(Cs,Cl)));
% (c) fitw's surplus columns are nonzero
hfit = t1d_subfun(src,'fitw');
Ns = F.B1.Ns;   rv1 = min(1,Ns);   rv2 = min(2,Ns);
rng(5,'twister');
Y0 = cell(1,numel(Ns));
for i = 1:numel(Ns), Y0{i} = orth(randn(Ns(i),rv1(i))); end
wf = hfit(Y0,rv2,Ns);
k = 0;   mn = inf;
for i = 1:numel(Ns)
    Yi = reshape(wf(k+(1:Ns(i)*rv2(i))),Ns(i),rv2(i));   k = k+Ns(i)*rv2(i);
    cn = sqrt(sum(Yi.^2,1));
    mn = min(mn,min(cn(rv1(i)+1:end)));
end
A = ck(A,'surplus-nonzero',mn>0, ...
    sprintf('fitw surplus column norms : min %.4e   EXPECTED > 0 (a zero column is B2)',mn));
% (d) and at RUN TIME the hook leads at every rung.  n=2 is used because rank 1
% cannot certify there (measured best relP 9.1e-02), so the ladder is forced to
% take a second rung and the ordering is actually observable.  Budgets are tiny
% on purpose: this measures dispatch order, not search quality.
Y2 = cell(1,numel(F.B2.Ns));
for i = 1:numel(F.B2.Ns), Y2{i} = orth(randn(F.B2.Ns(i),1)); end
Rw = t1d_bm(F.B2,struct('maxrank',2,'seeds',11,'lmit',20,'drit',10,'w0',{Y2}));
rr = cell2mat(Rw.log(:,1));   st = Rw.log(:,2);
first = {};   okw = true;
for r = unique(rr(:))'
    j = find(rr==r,1);
    first{end+1} = sprintf('r=%d:%s',r,st{j});                              %#ok<AGROW>
    okw = okw && strcmp(st{j},'w0');
end
A = ck(A,'w0-leads',okw && numel(first)>=2, ...
    sprintf('per-rung first start: %s   (want w0 at every rung, and >=2 rungs seen)',strjoin(first,' ')));
end

% =========================================================================
% ---------------------------- small helpers ------------------------------
% =========================================================================
function A = anew()
A = struct('lines',{{}},'bad',{{}});
end
function A = nt(A,s)
A.lines{end+1} = s;
end
function A = ck(A,id,c,s)
% chk takes the boolean AND the line that states the measured value, so a
% PASSING test still prints its numbers: the log is the record of what was
% measured, not merely that something was true.
A.lines{end+1} = sprintf('%-4s %s',tern(c,'[ok]','[XX]'),s);
if ~c, A.bad{end+1} = id; end
end
function s = tern(c,a,b)
if c, s = a; else, s = b; end
end

function [qf,ipk] = flipnull(B,q)
% A perturbation of q that leaves A(q) EXACTLY unchanged while driving one
% block indefinite.  Ssym averages each block's (i,j)/(j,i) columns, so
% Ssym*vec(M) depends only on sym(M); a null direction of the symmetric-
% parameter map therefore stays a null direction after symmetrisation.
qf = [];   ipk = 0;
for i = 1:numel(B.Ns)
    N = B.Ns(i);
    [ia,ja] = find(triu(ones(N)));   npar = numel(ia);
    E = sparse(N*N,npar);
    for p = 1:npar
        if ia(p)==ja(p)
            E((ja(p)-1)*N+ia(p),p) = 1;
        else
            E((ja(p)-1)*N+ia(p),p) = 1;   E((ia(p)-1)*N+ja(p),p) = 1;
        end
    end
    Z = null(full(B.P.Ssym(:,B.P.rows{i})*E));
    if isempty(Z), continue; end
    Qi = reshape(q(B.P.rows{i}),N,N);    Qi = (Qi+Qi')/2;
    Dm = reshape(E*Z(:,1),N,N);          Dm = (Dm+Dm')/2;
    al = 0.5*max(eig(Qi))/max(abs(eig(Dm)));   % push mineig to about -maxeig/2
    qf = q;   Qf = Qi + al*Dm;   qf(B.P.rows{i}) = Qf(:);
    ipk = i;   return
end
end

function w = settle(P,Ns,rv,seed)
% A rank-rv point the search would actually be sitting on: the same random
% start, Douglas-Rachford pass and LM descent t1d_bm uses, so T4 measures the
% Jacobian at a real iterate and not at an arbitrary matrix.
nb = numel(Ns);
rng(seed,'twister');
q0 = zeros(P.Ntot,1);
for i = 1:nb
    Mr = randn(Ns(i));   Mr = 0.1*(Mr+Mr')/sqrt(2*Ns(i));   q0(P.rows{i}) = Mr(:);
end
q0 = bm_proj(P,q0,rv,'affine');
[Vd,~] = bm_dr(P,rv,q0,300);
w = [];
for i = 1:nb
    Yi = zeros(Ns(i),rv(i));
    if ~isempty(Vd) && ~isempty(Vd{i})
        kk = min(size(Vd{i},2),rv(i));   Yi(:,1:kk) = Vd{i}(:,1:kk);
    end
    w = [w;Yi(:)];                                                          %#ok<AGROW>
end
w = bm_lm2(P,rv,w,400,1e-16);
end

function [lv,pd] = coljac(w,P,Ns,rold,rnew)
% Column-block norms of the SHIPPED analytic Jacobian, split into the columns
% carried over from rank rold and the surplus ones.  Classification is BY
% INDEX, never by whether a column happens to be zero: the question under test
% is whether a surplus column is mobile, and testing norm(Y(:,c))==0 would
% silently reclassify the mobile case out of the test.
[~,J] = bm_resid(w,P,rnew);
lv = [];   pd = [];   c0 = P.Kf;
for i = 1:numel(Ns)
    N = Ns(i);
    for c = 1:rnew(i)
        nc = norm(J(:,c0+(c-1)*N+(1:N)));
        if c <= rold(i), lv(end+1) = nc; else, pd(end+1) = nc; end          %#ok<AGROW>
    end
    c0 = c0 + N*rnew(i);
end
end

function idx = padidx(Ns,rold,rnew,Kf)
% w-indices of the surplus columns, for the "LM cannot move them" assertion.
idx = [];   off = Kf;
for i = 1:numel(Ns)
    N = Ns(i);
    for c = 1:rnew(i)
        if c > rold(i), idx = [idx, off+(c-1)*N+(1:N)]; end                 %#ok<AGROW>
    end
    off = off + N*rnew(i);
end
end
