function R = t1d_bm(B,o)                                                    % CC, 09/22/2026
% PROVENANCE.  scratchpad/reach1d/r1bm.m, renamed (B6), plus one input and two
% outputs added for the regression suite.  The search loop itself is untouched.
%   .w0      the seeding hook, dispatched EXACTLY as pielr_certify dispatches
%            it (w0 first at every rung, then the seeds, then 'prev'), using a
%            verbatim copy of its fitw.  T9 asserts that ordering here at run
%            time AND on pielr_certify's own source, and asserts the fitw copy
%            still matches, because the shipped hook lives in a 2-D-only
%            executive this suite cannot afford to run.
%   R.log    one row per attempt {rank, start name}: lets a test assert the
%            ORDER starts were tried in rather than trusting a comment.
%   R.V,R.w  the accepted point's orthonormal face and its factor, which T5
%            needs before it can ask whether a WIDER face still certifies.
% .lmit was already an option.  It is called out here only because the SHIPPED
% default is 400 and T6 is the test that exists because 400 is a hard cut, not
% a convergence test (defect B1).
% bm_setup/bm_dr/bm_lm2/bm_proj/bm_resid/bm_report/raw_data are the SHIPPED
% files, with byte-identity asserted by T0.
%                                                                 % CC, 09/22/2026
% r1bm(B,o) -- CURVE B: the SHIPPED Burer-Monteiro discovery, in 1-D.
%
% The loop below is pielr_certify's disc_bm with the 2-D calls swapped for
% their 1-D equivalents: rank ladder r = 1..maxrank, per rank a random start
% per seed (initial magnitude picked by the seed INDEX from
% [1e-1 1 1e1 1e2 1e4 1e6]) taken through 300 Douglas-Rachford iterations and
% then bm_lm2, plus a zero-padded warm start from the best point at the
% previous rank; accept as soon as EITHER the BM point itself or the
% face-restricted re-solve passes the operator gate.  bm_setup, bm_proj, bm_dr,
% bm_resid, bm_report and bm_lm2 are the shipped files unmodified (verified
% bit-identical by diff against PIETOOLS_demos/lowrank_2d_stability/private);
% restrict_solve_1d differs from the shipped restrict_solve in one line, the
% gate call.  So this measures the shipped search, not a re-implementation.
%
% COST IS REPORTED MACHINE-INDEPENDENTLY FIRST: attempts, Levenberg-Marquardt
% iterations, Douglas-Rachford iterations, gate evaluations.  Wall clock is
% recorded but is a weaker number.
if nargin<2, o = struct(); end
dflt = struct('maxrank',8,'seeds',[11 22 33 44],'drit',300,'lmit',400, ...
              'lmtol',1e-16,'verb',0,'timecap',inf,'w0',[]);                % CC, 09/22/2026
fn = fieldnames(dflt);
for i=1:numel(fn), if ~isfield(o,fn{i})||isempty(o.(fn{i})), o.(fn{i})=dflt.(fn{i}); end, end
P = B.P;   Ns = B.Ns;   nb = numel(Ns);
R = struct('cert',false,'route','','r',[],'rank',[],'relP',NaN,'relD',NaN, ...
           'mineig',[],'attempts',0,'lm_iters',0,'dr_iters',0,'gates',0, ...
           'faces',0,'sdp_faces',0,'ranks_tried',[],'best_relP',inf,'best_rawrel',inf, ...
           'best_r',[],'t',NaN,'seed_hit',NaN,'start_hit','','err','', ...
           'per_rank',[],'q',[],'log',{{}},'V',{{}},'w',[]);                % CC, 09/22/2026
% R.V: the orthonormal FACE of the accepted point, and R.w its factor.  T5
% needs a face that is known to certify before it can test whether a WIDER
% face still does (defect B3), and manufacturing one from the interior-point
% solution would test a different object than the search produces.  % CC, 09/22/2026
% R.q: the Gram vector of the ACCEPTED point, normalised-b units, so the
% acceptance can be re-run afterwards with the Gram-free Galerkin
% confirmations switched on (gate1d's do_gal).  Purely an extra output -- it
% consumes no randomness and does not touch the search.
T0 = tic;   prevw = [];   prevrv = [];
PR = [];
for r = 1:o.maxrank
    rv = min(r,Ns);
    R.ranks_tried(end+1) = r;                                        %#ok<AGROW>
    bestraw = inf;  bestw = [];  bestrelP_r = inf;
    % CC, 09/22/2026: dispatch copied from pielr_certify's disc_bm -- a
    % supplied .w0 leads at EVERY rung.  It has to lead: on this family seed 11
    % certifies on the first attempt at low lam, so an arm appended after the
    % seeds would never execute and the option would silently measure nothing
    % (defect B7).
    nw0 = double(~isempty(o.w0));                                           % CC, 09/22/2026
    nst = nw0 + numel(o.seeds) + ~isempty(prevw);                           % CC, 09/22/2026
    for si = 1:nst
        if nw0 && si==1                                                     % CC, 09/22/2026
            w = fitw(o.w0,rv,Ns);   sname = 'w0';   shit = NaN;             % CC, 09/22/2026
        elseif si <= nw0 + numel(o.seeds)                                   % CC, 09/22/2026
            sj = si - nw0;                                                  % CC, 09/22/2026
            rng(o.seeds(sj),'twister');                                     % CC, 09/22/2026
            scl = [1e-1 1 1e1 1e2 1e4 1e6];  sc = scl(mod(sj-1,numel(scl))+1); % CC, 09/22/2026
            q0 = zeros(P.Ntot,1);
            for i=1:nb
                Mr = randn(Ns(i));  Mr = sc*(Mr+Mr')/sqrt(2*Ns(i));
                q0(P.rows{i}) = Mr(:);
            end
            q0 = bm_proj(P,q0,rv,'affine');
            [Vd,~] = bm_dr(P,rv,q0,o.drit);
            R.dr_iters = R.dr_iters + o.drit;
            w = [];
            for i=1:nb
                Yi = zeros(Ns(i),rv(i));
                if ~isempty(Vd)&&~isempty(Vd{i})
                    kk = min(size(Vd{i},2),rv(i));  Yi(:,1:kk) = Vd{i}(:,1:kk);
                end
                w = [w;Yi(:)];                                       %#ok<AGROW>
            end
            sname = sprintf('seed%d',o.seeds(sj));   shit = o.seeds(sj);    % CC, 09/22/2026
        else
            w = padw(prevw,prevrv,rv,Ns);   sname = 'prev';   shit = NaN;
        end
        % the order starts were tried in, for T9's dispatch assertion
        R.log(end+1,:) = {r,sname};                                         % CC, 09/22/2026
        [w,~,itlm] = bm_lm2(P,rv,w,o.lmit,o.lmtol);
        R.attempts = R.attempts + 1;   R.lm_iters = R.lm_iters + itlm;
        R1 = bm_report(w,P,rv);
        if R1.raw_rel < bestraw, bestraw = R1.raw_rel;  bestw = w; end
        if R1.raw_rel < R.best_rawrel, R.best_rawrel = R1.raw_rel; end
        [~,~,qk] = bm_resid(w,P,rv);
        Gk = gate1d(B.prog,B.H,P,qk);   R.gates = R.gates + 1;
        if Gk.relP < R.best_relP, R.best_relP = Gk.relP; R.best_r = rv; end
        if Gk.relP < bestrelP_r, bestrelP_r = Gk.relP; end
        % certify the SUBSPACE the search found -- the find/certify split
        Vc = cell(1,nb);  k = 0;
        for i=1:nb
            Yi = reshape(w(k+(1:Ns(i)*rv(i))),Ns(i),rv(i));  k = k+Ns(i)*rv(i);
            Vi = orth(Yi);
            if isempty(Vi), Vi = zeros(Ns(i),1); Vi(1)=1; end
            Vc{i} = Vi;
        end
        % Face certification, BOTH arms.  The shipped default solves the
        % restricted system by backslash on a pivoted-QR row selection, which
        % is the intended DETERMINED solve only while np = sum r(r+1)/2 does
        % not exceed rM = rank(M) <= m.  Measured here (smoke_F2): at 'heavy'
        % m = 56 and np overtakes rM at r = 5, after which backslash returns a
        % basic solution and the arm rejects even the IPM certificate's OWN
        % face (r = [10 17 10], relP 5.9 with clip 2.5) -- a false negative.
        % The SeDuMi arm accepts that face (relP 5.0e-7).  Both are sound
        % one-way, so trying the second after the first only makes BM's
        % measured reach more generous, never falsely larger.
        Gf = [];  frel = NaN;  fsdp = false;  qf = [];
        try
            [Gf,qf] = restrict_solve_1d(B.prog,B.H,P,B.Atf,B.bf,Vc,false);
            R.faces = R.faces + 1;   frel = Gf.relP;
            if Gf.relP < R.best_relP, R.best_relP = Gf.relP; R.best_r = rv; end
            if Gf.relP < bestrelP_r, bestrelP_r = Gf.relP; end
        catch ME
            R.err = ME.message;
        end
        if isempty(Gf) || ~Gf.cert
            try
                [G2,q2] = restrict_solve_1d(B.prog,B.H,P,B.Atf,B.bf,Vc,true);
                R.faces = R.faces + 1;   R.sdp_faces = R.sdp_faces + 1;
                if G2.relP < R.best_relP, R.best_relP = G2.relP; R.best_r = rv; end
                if G2.relP < bestrelP_r, bestrelP_r = G2.relP; end
                if isempty(Gf) || G2.relP < Gf.relP
                    Gf = G2;  frel = G2.relP;  fsdp = true;  qf = q2;
                end
            catch ME
                R.err = ME.message;
            end
        end
        if o.verb
            fprintf('    r=%-2d %-8s raw_rel %9.3e  BM relP %9.3e  face relP %9.3e\n', ...
                r,sname,R1.raw_rel,Gk.relP,frel);
        end
        if ~isempty(Gf) && Gf.cert
            R.cert=true; R.r=rv; R.rank=Gf.rank;
            if fsdp, R.route='face-sdp'; else, R.route='face-ls'; end
            R.relP=Gf.relP; R.relD=Gf.rel; R.mineig=Gf.mineig;  R.q=qf;
            R.V=Vc;  R.w=w;                                                 % CC, 09/22/2026
            R.seed_hit=shit; R.start_hit=sname;  R.t=toc(T0);
            R.per_rank = PR;  return
        elseif Gk.cert
            R.cert=true; R.route='bm'; R.r=rv; R.rank=Gk.rank;
            R.relP=Gk.relP; R.relD=Gk.rel; R.mineig=Gk.mineig;  R.q=qk;
            R.V=Vc;  R.w=w;                                                 % CC, 09/22/2026
            R.seed_hit=shit; R.start_hit=sname;  R.t=toc(T0);
            R.per_rank = PR;  return
        end
        if toc(T0) > o.timecap, R.err='timecap'; R.t=toc(T0); R.per_rank=PR; return, end
    end
    PR = [PR; r bestraw bestrelP_r R.attempts R.lm_iters];           %#ok<AGROW>
    prevw = bestw;  prevrv = rv;
end
R.t = toc(T0);   R.per_rank = PR;
end

function w = fitw(Y0,rv,Ns)                                                 % CC, 09/22/2026
% VERBATIM copy of pielr_certify's fitw subfunction (that file is 2-D-only, and
% MATLAB gives no way to call a subfunction from outside).  T9 asserts this
% copy still matches the shipped one after the change markers are stripped; if
% it ever does not, the suite has stopped testing the shipped fill rule.
% Fit a supplied per-block factor {Y_1..Y_B} to the current rank profile rv.
% Keep the leading columns; fill any surplus with SMALL RANDOM values scaled
% to the factor's own entries.  Never zeros: a zero column of Y has an
% identically zero Jacobian block (J_i = 2*W*Ssym(:,rows_i)*kron(Y_i,I), whose
% column block for column c is Y_i(:,c) (x) I), so LM can never move it and the
% extra rank would be inert -- the same fixed point the 'prev' pad sits on.
w = [];
for i = 1:numel(Ns)
    Yi = zeros(Ns(i),rv(i));
    Y  = Y0{i};
    kk = min(size(Y,2),rv(i));
    if kk>0, Yi(:,1:kk) = Y(:,1:kk); end
    if rv(i) > kk
        s = norm(Y,'fro')/max(sqrt(numel(Y)),1);       % rms entry of the seed
        if ~isfinite(s) || s<=0, s = 1; end
        Yi(:,kk+1:end) = 1e-3*s*randn(Ns(i),rv(i)-kk);
    end
    w = [w;Yi(:)]; %#ok<AGROW>
end
end

function w = padw(wold,rold,rnew,Ns)
w=[]; k=0;
for i=1:numel(Ns)
    Yo = reshape(wold(k+(1:Ns(i)*rold(i))),Ns(i),rold(i)); k=k+Ns(i)*rold(i);
    Yn = zeros(Ns(i),rnew(i)); c=min(rold(i),rnew(i)); Yn(:,1:c)=Yo(:,1:c);
    w=[w;Yn(:)];                                                     %#ok<AGROW>
end
end
