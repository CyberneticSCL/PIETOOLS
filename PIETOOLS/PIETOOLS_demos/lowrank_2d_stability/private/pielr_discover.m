function [V,R,q,why,notes] = pielr_discover(prog,H,D,P,A,opts,vb)          % CC, 09/23/2026
% PIELR_DISCOVER  Burer-Monteiro face discovery and face certification, for
% any adapted executive and either dimension.
%
% Same pipeline as pielr_certify's local disc_bm and, apart from the gate
% (see pielr_opcheck), the same arithmetic: a rank ladder; per rung a
% Douglas-Rachford start per seed plus the previous rung's point; LM on the
% whitened equality residual; then gate the BM point itself, and otherwise
% hand orth(Y_i) to restrict_solve as a face.
%
% WHAT IS NEW HERE IS THAT NOTHING IS THROWN AWAY.  Three diagnostics were
% being computed on every attempt and discarded at the call site:
%   * bm_lm2 returns WHY it stopped ('maxit' | 'tol' | 'stagnant' |
%     'damping').  Every existing call site takes one output.  Without it a
%     run that hit its iteration cap is indistinguishable from one that
%     converged, which is exactly the question the package's own open
%     questions about budget-vs-rank turn on.
%   * restrict_solve returns info.np, info.rM and info.rel_eq_full.  Every
%     existing call site takes two outputs.  rM in particular is recorded as
%     "never measured on the 2-D program" while being computed on every face.
%   * the gate returns the maxima its ratio is built from.  Only the ratio was
%     kept, so no log can say whether two residuals differ in the numerator or
%     in the denominator.
% All three go into the per-attempt table WHY, which is what makes one run
% comparable with another.
%
% NOTE ON THE SEEDS, worth reading before comparing attempts: the seed INDEX
% selects the initial magnitude from [1e-1 1 1e1 1e2 1e4 1e6], so seeds
% [11 22 33] are a designed sweep over three scales, not three independent
% samples of one distribution.
%
% OUTPUT  V,R,q  accepted face, gate result and Gram vector ([] if none)
%         why    struct array, one row per attempt
%         notes  cellstr of fine print

notes = {};  V = [];  R = [];  q = [];
Ns = P.Ns;   B = numel(Ns);
why = struct('rank',{},'start',{},'raw_rel',{},'bm_rel',{},'bm_rel_d',{}, ...
             'lm_exit',{},'lm_iters',{},'face_rel',{},'np',{},'rM',{}, ...
             'rel_eq_full',{},'clip',{},'verdict',{});

% rank ladder: opts.rank first (if given), then upward to maxrank
rvl = {};
if ~isempty(opts.rank)
    rv0 = opts.rank(:)';
    if isscalar(rv0), rv0 = rv0*ones(1,B); end
    assert(numel(rv0)==B,'opts.rank must be scalar or one entry per block');
    rvl{end+1} = min(rv0,Ns);
    rnext = max(rv0)+1;
else
    rnext = 1;
end
for r = rnext:opts.maxrank, rvl{end+1} = min(r,Ns); end %#ok<AGROW>

if vb
    fprintf('  %-11s | %-8s | %-10s | %-10s | %-8s | %-10s | %s\n', ...
        'rank','start','BM raw','BM rel','LM exit','face rel','verdict');
end

prevw = [];  prevrv = [];
for ri = 1:numel(rvl)
    rv = rvl{ri};
    bestraw = inf;  bestw = [];
    nw0 = double(~isempty(opts.w0));
    nst = nw0 + numel(opts.seeds) + ~isempty(prevw);
    for si = 1:nst
        if nw0 && si==1
            w = pielr_fitw(opts.w0,rv,Ns,zeros(P.Kf,1));   sname = 'w0';
        elseif si <= nw0 + numel(opts.seeds)
            sj = si - nw0;
            rng(opts.seeds(sj),'twister');
            scl = [1e-1 1 1e1 1e2 1e4 1e6];  sc = scl(mod(sj-1,numel(scl))+1);
            q0 = zeros(P.Ntot,1);
            for i = 1:B
                Mr = randn(Ns(i));  Mr = sc*(Mr+Mr')/sqrt(2*Ns(i));
                q0(P.rows{i}) = Mr(:);
            end
            q0 = bm_proj(P,q0,rv,'affine');
            [Vd,~] = bm_dr(P,rv,q0,300);
            % the free coordinates come from the affine projection, which has
            % already put a consistent value there; disc_bm omitted them
            % entirely, which is correct only when Kf == 0
            w = q0(P.free);                                                 % CC, 09/23/2026
            for i = 1:B
                Yi = zeros(Ns(i),rv(i));
                if ~isempty(Vd) && ~isempty(Vd{i})
                    kk = min(size(Vd{i},2),rv(i));  Yi(:,1:kk) = Vd{i}(:,1:kk);
                end
                w = [w;Yi(:)]; %#ok<AGROW>
            end
            sname = sprintf('seed%d',opts.seeds(sj));
        else
            % Warm start from the previous rung.  pielr_fitw, not the zero pad:
            % a zero column of Y makes its whole Jacobian column block
            % 2*W*Ssym(:,rows)*kron(Y(:,c),I) identically zero, so the LM step
            % has no component there and the added rank can never move.  The
            % package's own T4 pins that as defect B2 and fitw is its tested
            % repair; the zero pad is left in pielr_certify, which remains the
            % measurement baseline.
            Y0 = cell(1,B);  k = P.Kf;      % w starts with the free block
            for i = 1:B
                Y0{i} = reshape(prevw(k+(1:Ns(i)*prevrv(i))),Ns(i),prevrv(i));
                k = k + Ns(i)*prevrv(i);
            end
            w = pielr_fitw(Y0,rv,Ns,prevw(1:P.Kf));   sname = 'prev';
        end

        [w,~,itlm,lmexit] = bm_lm2(P,rv,w,opts.lmit,opts.lmtol);
        R1 = bm_report(w,P,rv);
        if R1.raw_rel < bestraw, bestraw = R1.raw_rel;  bestw = w; end
        [~,~,qk] = bm_resid(w,P,rv);
        Rk = pielr_opcheck(prog,H,P,qk,A);

        % certify the SUBSPACE BM found -- the point of the find/certify split
        Vc = cell(1,B);
        for i = 1:B
            Yi = reshape(pielr_wblk(w,Ns,rv,i,P.Kf),Ns(i),rv(i));
            Vi = orth(Yi);
            if isempty(Vi)                 % Y_i == 0: keep a 1-dim face, which
                Vi = zeros(Ns(i),1); Vi(1) = 1;   % still contains X_i = 0
            end
            Vc{i} = Vi;
        end
        Rf = [];  qf = [];  frel = NaN;  inf_ = struct('np',NaN,'rM',NaN, ...
                                          'rel_eq_full',NaN,'clip',NaN);
        try
            [Rf,qf,inf_] = restrict_solve(prog,H,P,D.At,D.b,Vc,false,A);
            frel = Rf.rel;
        catch ME
            if vb, fprintf('      (restrict_solve: %s)\n',ME.message(1:min(70,end))); end
        end
        okf = ~isempty(Rf) && Rf.ok;
        okb = Rk.ok;
        vd  = 'fail';
        if okf, vd = 'FACE OK'; elseif okb, vd = 'BM OK'; end

        why(end+1) = struct('rank',rv,'start',sname,'raw_rel',R1.raw_rel, ...
            'bm_rel',Rk.rel,'bm_rel_d',Rk.rel_d,'lm_exit',lmexit, ...
            'lm_iters',itlm,'face_rel',frel,'np',getf(inf_,'np'), ...
            'rM',getf(inf_,'rM'),'rel_eq_full',getf(inf_,'rel_eq_full'), ...
            'clip',getf(inf_,'clip'),'verdict',vd); %#ok<AGROW>
        if vb
            fprintf('  %-11s | %-8s | %-10.4g | %-10.4g | %-8s | %-10.4g | %s\n', ...
                mat2str(rv),sname,R1.raw_rel,Rk.rel,lmexit,frel,vd);
        end

        if okf
            V = Vc;  R = Rf;  q = qf;  return
        elseif okb
            V = Vc;  R = Rk;  q = qk;
            notes{end+1} = ['accepted the BM point directly: the operator gate ' ...
                'passed without a face re-solve, so this point satisfies the ' ...
                'equalities only to the LM residual and the face machinery did ' ...
                'no work'];   %#ok<AGROW>
            return
        end
    end
    prevw = bestw;  prevrv = rv;
end
end

% =========================================================================
function v = getf(s,f)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = NaN; end
end
