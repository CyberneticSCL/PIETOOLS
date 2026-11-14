function D = heatNd_sdp(prog,file,meta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% D = HEATND_SDP(PROG,FILE,META) the SeDuMi-form SDP
%
%       min c'x   s.t.   At'x = b,   x in K = R^Kf x S^Ks(1)_+ x ...,
%
% that sossolve would hand the solver for the UNSOLVED program PROG, with
% nothing shadowed, optionally saved (-v7.3) with metadata for solver
% benchmarking. Replicates sossolve.m (SOSTOOLS400, working tree of
% 09/27/2026): the expression concatenation, 'processvars' (a local
% function of sossolve, so not callable; 'poly' -> spantiblkdiag, 'sos' ->
% spblkdiag), the objective c = RR'*c0, and the b normalization b/||b||
% (MMP 09/12/2026; x_orig = x*bscl). Same construction as the api scout's
% probe (hapi_sdp, not in the repository), whose re-solve reproduced
% lpisolve's Mosek status and residual.
%
% ONLY FOR PROGRAMS WITHOUT 'ineq' EXPRESSIONS: sossolve adds their slack
% Gram blocks during the solve (addextrasosvar), so a pre-solve dump would
% be incomplete; this errors instead. HEATND_LPI makes 'eq' rows only.
% Not replicated (errors never, since not requested): options.simplify /
% frlib reductions, which sossolve applies only when asked.
%
% INPUT
% - prog: unsolved program (HEATND_LPI), or an SDP struct D already built
%         by this function (e.g. from a k-affine family), to save only;
%         D.keep (row set, HEATND_LINDEP) and D.ref (reference verdict:
%         st, rows 'keep'|'all', tight, rel_b or cert_viol (cert_rel in
%         dumps written before the final reviews, 09/27/2026), t_mosek,
%         mosek_threads, iter, when, source) are saved in info when present,
%         so HEATND_SOLVE(FILE) repeats the route that produced D.ref;
% - file: (optional) .mat path; saves At, b, c, K, bscl, info (v7.3);
% - meta: (optional) HEATND_LPI meta, stored in info (base removed); its
%         kappa (the only rate the SDP sees) is copied to info.kappa.
% OUTPUT struct D: At (nx x m sparse), b (m x 1), c (nx x 1), K (f,l,q,s),
%   bscl, RR (cone order -> decvartable order: RRx = RR*x*bscl), m, nx,
%   nnz, zero_rows (rows with no variable), deficient (a zero row with b ~=
%   0: sossolve's presolve declares the program infeasible before any
%   solver runs, sossolve.m "trivially infeasible"; a DEGREE deficiency,
%   0 = b_i ~= 0, not a property of k), info (what was saved).
%
% Cost: one horizontal concatenation of the expression blocks and one
% RR'*Atf, O(nnz); saving is O(nnz) on disk (v7.3 compresses).
%
% Initial coding MMP, 09/27/2026
% MMP, 09/27/2026 (review fixes): save an SDP struct as given, with its
%   row set, reference verdict/route and kappa, so a dump is reproducible.
% MMP, 09/27/2026 (final reviews): header only - info.ref names cert_viol.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    file = '';      end
if nargin<3,    meta = [];      end
if isstruct(prog) && isfield(prog,'At')        % an SDP struct: save only
    D = prog;
    if ~isempty(file),  D = save_dump(D,file,meta);     end
    return
end
if ~isempty(prog.solinfo.x),    error('heatNd_sdp:solved','Program is already solved.'),  end
typ = prog.expr.type(1:prog.expr.num);
bad = ~strcmp(typ,'eq');
if any(bad)
    error('heatNd_sdp:ineq',['Expression %d is ''%s''; sossolve adds slack blocks during ' ...
          'the solve, so a pre-solve dump would be incomplete.'],find(bad,1),typ{find(bad,1)});
end
if isfield(prog,'extravar') && isfield(prog.extravar,'num') && prog.extravar.num>0
    error('heatNd_sdp:extravar','Program already carries slack variables.');
end
Atc = prog.expr.At(1:prog.expr.num);    bc = prog.expr.b(1:prog.expr.num);
Atf = [Atc{:}];     bf = vertcat(bc{:});                   % one concatenation
% processvars (sossolve.m): free decision variables first, then per variable.
RR = speye(prog.var.idx{1}-1);  Kf = prog.var.idx{1}-1;     Ks = [];
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly',    RR = spantiblkdiag(RR,speye(sz));   Kf = Kf+sz;
        case 'sos',     RR = spblkdiag(RR,speye(sz));       Ks(end+1) = round(sqrt(sz)); %#ok<AGROW>
        otherwise,      error('heatNd_sdp:vartype','Unexpected variable type ''%s''.',prog.var.type{i});
    end
end
At = RR'*Atf;
bscl = norm(bf);    if bscl==0 || ~isfinite(bscl),  bscl = 1;   end
b = bf/bscl;
c0 = sparse(size(Atf,1),1);
if ~isempty(prog.objective),    c0(1:numel(prog.objective)) = prog.objective;   end
c = RR'*c0;
K = struct('f',Kf,'l',0,'q',[],'s',Ks);
ptol = 1e-12;                                   % sossolve: max(pars.tol*1e-3,1e-12)
zr = full(max(abs(At),[],1))'<=ptol;            % sossolve's test; transient 16 B/nnz
D = struct('At',At,'b',b,'c',c,'K',K,'bscl',bscl,'RR',RR,'m',size(At,2), ...
           'nx',size(At,1),'nnz',nnz(At),'zero_rows',nnz(zr), ...
           'deficient',full(any(zr & abs(b)>ptol*max(abs(b)))));
if ~isempty(file),  D = save_dump(D,file,meta);     end
end


function D = save_dump(D,file,meta)
% Save D (v7.3) with info: shape, deficiency, meta, kappa, and - when D
% carries them - the independent row set and the reference verdict/route.
info = struct('what',['Jagt & Peet arXiv:2508.14840v4 Cor. 35, heatNd benchmark: ' ...
    'SeDuMi form min c''x s.t. At''*x = b, x in K (K.f free, then K.s PSD blocks, ' ...
    'column-major vec); b normalized, original x = x*bscl. The SDP depends on ' ...
    'kappa = r + k only.'], ...
    'm',size(D.At,2),'nx',size(D.At,1),'nnz',nnz(D.At),'deficient',D.deficient, ...
    'created',char(datetime('now')),'kappa',NaN);   % sizes from At: D may be a k-affine combination
if ~isempty(meta)
    mt = meta;  if isfield(mt,'base'),  mt = rmfield(mt,'base');   end
    info.meta = mt;
    if isfield(mt,'kappa') && ~isempty(mt.kappa),   info.kappa = mt.kappa;  end
end
if isfield(D,'kappa') && ~isempty(D.kappa),     info.kappa = D.kappa;   end
if isfield(D,'keep') && ~isempty(D.keep),       info.keep = D.keep(:);  end
if isfield(D,'ref') && ~isempty(D.ref),         info.ref = D.ref;       end
At = D.At;  b = D.b;    c = D.c;    K = D.K;    bscl = D.bscl;  %#ok<NASGU> saved below
save(file,'At','b','c','K','bscl','info','-v7.3');
D.info = info;
end
