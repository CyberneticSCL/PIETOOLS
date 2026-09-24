function ref = pielr_ipm_ref(prog,H,D,A,vb,solver,rb)                       % CC, 09/23/2026
% PIELR_IPM_REF  Solve the SAME assembled program with the interior-point
% method and judge it with the SAME gate, to give the low-rank result a
% reference measured in this session.
%
% WHY THE SAME PROGRAM AND THE SAME GATE.  Comparing a low-rank certificate
% against a number quoted in a library comment compares two things that were
% produced by different code at different settings; and comparing against a
% solver's own feasibility flag compares a residual to a status word.  Here
% the interior-point solver is handed the identical prog that the low-rank
% search was handed, its solution is pushed through pielr_opcheck exactly as
% the candidate is, and the two are then commensurable by construction.
%
% For an LPI with an objective the solver minimises it, so ref.gam is the
% quantity the low-rank bisection is trying to approach from above.
%
% OUTPUT  ref.ok    the gate's verdict on the IPM point
%         ref.rel .rel_d .maxRes .maxNrm .mineig .rank   as pielr_opcheck
%         ref.gam  objective value (NaN when the LPI has none)
%         ref.t    seconds
%         ref.err  error message if the solve failed, '' otherwise

ref = struct('ok',false,'rel',NaN,'rel_d',NaN,'maxRes',NaN,'maxNrm',NaN, ...
             'solver','','rel_obj',NaN,'fixed',false, ...
             'mineig',NaN,'rank',NaN,'gam',NaN,'t',NaN,'err','');
% THE REFERENCE SOLVER IS AN ARGUMENT, and it is pinned.  sossolve takes the
% first solver on the path and mosekopt leads its list, so an unset
% sos_opts.solver silently selects Mosek -- verified present on this machine
% at C:\Program Files\Mosek\11.0.  A reference that does not say which solver
% produced it is not a reference.
%
% 'best' (the default) SOLVES WITH EVERY AVAILABLE SOLVER AND KEEPS THE BEST
% RESULT, because neither dominates and which one wins is program-dependent.
% MEASURED on 2-D stability, identical programs, eppos 1e-2, gate quantity
% max|Res|/max|Pop|:
%   deg frac   SeDuMi                Mosek                 better
%    3  0.10   1.8996e-06            7.1431e-05            SeDuMi  27x
%    3  0.25   5.1666e-06            1.5419e+00            SeDuMi  3e5x
%    3  0.50   6.2832e-05            5.2404e-03            SeDuMi  83x
%    4  0.10   7.0682e-07  (ok)      4.7985e-09  (ok)      Mosek  147x
%    4  0.25   1.4780e-06  (fail)    2.8121e-08  (ok)      Mosek   53x
%    4  0.50   2.6697e-06            5.0261e+00            SeDuMi  2e6x
% Mosek wins by two orders at degree 4 on the easier operating points and
% CERTIFIES at frac 0.25 where SeDuMi does not; SeDuMi wins by up to six
% orders everywhere at degree 3 and at degree 4 frac 0.50.  Picking either a
% priori mismeasures half the table, and since the score rel/ipm_rel means
% "how close to the BEST AVAILABLE answer", the reference has to actually be
% the best available.  Mosek is also 4-6x faster, so the extra solve is cheap
% relative to the low-rank arm it is scoring.
if nargin<6 || isempty(solver), solver = 'best'; end                        % CC, 09/23/2026
if nargin<7, rb = []; end                                                   % CC, 09/24/2026
if strcmpi(solver,'best')
    cand = {};
    if ~isempty(which('sedumi')),   cand{end+1} = 'sedumi';   end
    if ~isempty(which('mosekopt')), cand{end+1} = 'mosek';    end
    if isempty(cand), ref.err = 'no interior-point solver on the path'; return, end
    best = [];
    for i = 1:numel(cand)
        % rb is threaded through, so each candidate is ranked on the SAME
        % quantity the score will use -- the fixed-gamma residual.  Ranking
        % on the optimising residual and then re-solving the winner would
        % pick a solver by a number nothing downstream reads. % CC, 09/24/2026
        ri = pielr_ipm_ref(prog,H,D,A,vb,cand{i},rb);
        if isempty(ri.err) && (isempty(best) || ri.rel < best.rel), best = ri; end
    end
    if isempty(best)
        ref.err = 'every available solver failed';  return
    end
    ref = best;   return
end
o.solver = solver;   o.simplify = false;
t0 = tic;
try
    if vb
        ps = lpisolve(prog,o);
    else
        evalc('ps = lpisolve(prog,o);');
    end
catch ME
    ref.err = ME.message;   ref.t = toc(t0);   return
end
ref.t = toc(t0);   ref.solver = solver;

% RRx, NOT x.  solinfo.x is SeDuMi's primal for the cone the solver actually
% saw: sossolve turns each 'ineq' expression into a slack cone variable, so on
% the 1-D l2gain program x has 498 entries against the program's 497 decision
% variables, and its ordering is the solver's, not the program's.  The rows of
% expr.At -- which is what bm_setup, the layout and the gate all index against
% -- are in DECVARTABLE order, and solinfo.RRx is the decision vector in that
% order.  MEASURED, using x: the executive's own solution violates its own
% equality rows by 1.83e+01 and x(jc) reads 2.562 for a program whose reported
% gamma is 0.516333.  The two coincide only when the program has no free
% variables and no inequality, which is exactly the stability case that had
% been the only thing tested.
x = full(ps.solinfo.RRx(:));
% sossolve APPENDS a slack decision variable per 'ineq' expression, to both the
% cone and decvartable, so the solved RRx is longer than the pre-solve layout
% by exactly that many entries and the program's own coordinates are the
% LEADING Ntot.  MEASURED on the 1-D l2gain program: 498 against 497, one
% inequality (gam >= 0), and truncating to 497 makes the executive's own
% solution pass the gate at rel = 1.01e-08.
nslack = numel(x) - D.L.Ntot;
if nslack < 0 || nslack > size(D.Gineq,2)
    ref.err = sprintf(['RRx has %d entries, layout says %d, and there are ' ...
        '%d inequality row(s) -- unexplained difference'], ...
        numel(x),D.L.Ntot,size(D.Gineq,2));
    return
end
x = x(1:D.L.Ntot);
% the gate wants NORMALISED-b units; build the bm_setup package only for the
% scaling and the block rows, so pre=0 (no whitener is needed to gate a point)
P = bm_setup(D.At,D.b,D.L.N,D.L.Kf,0,D.L);
R = pielr_opcheck(prog,H,P,x/P.nb0,A);
ref.ok = R.ok;   ref.rel = R.rel;   ref.rel_d = R.rel_d;
ref.maxRes = R.maxRes;  ref.maxNrm = R.maxNrm;
ref.mineig = R.mineig;  ref.rank = R.rank;
jc = find(D.c);
if numel(jc)==1, ref.gam = x(jc); end

% ---- MATCH THE COMPUTATION THE CANDIDATE ACTUALLY PERFORMS --------------
% (CC, 09/24/2026) When the LPI carries an objective, the two arms were not
% solving the same kind of problem: pielr_bisect_obj never optimises -- it
% PINS gamma and tests feasibility -- while this reference minimised.  With
% score = rel/ipm_rel the reference residual sets the acceptance threshold,
% so the mismatch is not cosmetic: an optimising solve reaches a worse
% equality residual (the iterates trade feasibility against optimality;
% with c = 0 the dual condition is satisfiable for free and the solver
% spends everything on primal feasibility), which INFLATES the threshold and
% admits candidate points that a like-for-like reference would reject.
%
% Given a rebuild handle, re-solve at the optimum as a pure feasibility
% program and report THAT residual, keeping gamma from the optimising solve.
% Without one the behaviour is exactly as before, so every existing caller
% is unaffected.
if nargin>=7 && ~isempty(rb) && isfinite(ref.gam)
    try
        [pg,Hg] = rb(ref.gam);
        Dg = pielr_rawdata(pg);
        rg = pielr_ipm_ref(pg,Hg,Dg,A,vb,solver);      % no rb: no recursion
        if isempty(rg.err) && isfinite(rg.rel)
            ref.rel_obj = ref.rel;        % keep both; the ratio is the finding
            ref.rel = rg.rel;   ref.ok = rg.ok;
            ref.maxRes = rg.maxRes;  ref.maxNrm = rg.maxNrm;
            ref.rel_d  = rg.rel_d;
            ref.mineig = rg.mineig;  ref.rank = rg.rank;
            ref.fixed  = true;
        end
    catch ME
        ref.err = ['fixed-gamma reference failed: ' ME.message];
    end
end
end
