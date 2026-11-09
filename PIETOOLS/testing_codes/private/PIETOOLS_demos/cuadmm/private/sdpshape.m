function S = sdpshape(prog)
% sdpshape(prog) -- the SDP a PIETOOLS program hands the solver.
%
% Replicates sossolve's own processvars/objective logic (sossolve.m:196-223)
% rather than calling it, because processvars is a local function there.  What
% is reproduced, and why each matters for a first-order solver:
%   K.f   free block     -- cuADMM's 'u' type; nonzero exactly when an executive
%                           declares a decision variable such as gam
%   K.s   PSD blocks     -- cuADMM does a DENSE eigendecomposition per block per
%                           iteration, so the block SIZE DISTRIBUTION, not m, is
%                           the per-iteration cost driver: sum(N^3)
%   c     objective      -- zero for the stability/well-posedness executives and
%                           a unit vector on gam for every Hinf/H2 one.  With
%                           c=0 the duality gap is degenerate and says nothing.
%   nnz   sparsity       -- cuADMM's memory is linear in nnz(At)

Atf = [];  bf = [];
for i = 1:prog.expr.num
    Atf = [Atf, prog.expr.At{i}];   %#ok<AGROW>
    bf  = [bf;  prog.expr.b{i}];    %#ok<AGROW>
end

% --- cone, in sossolve's order: decvars, then var (poly->free, sos->PSD),
%     then extravar (always PSD).
Kf = prog.var.idx{1}-1;
Ks = [];
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly', Kf = Kf + sz;
        case 'sos',  Ks(end+1) = sqrt(sz);   %#ok<AGROW>
    end
end
for i = 1:prog.extravar.num
    Ks(end+1) = sqrt(prog.extravar.idx{i+1}-prog.extravar.idx{i});  %#ok<AGROW>
end

% --- objective.  sossolve puts sos.objective on the leading block of c.
nvar = size(Atf,1);
cnz  = 0;
if isfield(prog,'objective') && ~isempty(prog.objective)
    cnz = nnz(prog.objective);
end

S.m        = numel(bf);
S.nvar     = nvar;
S.Kf       = Kf;
S.Ks       = Ks(:)';
S.nblk     = numel(Ks);
S.Nmax     = max([Ks 0]);
S.Nsum     = sum(Ks);
S.eigcost  = sum(double(Ks).^3);     % dense eig work per ADMM iteration
S.svec_len = Kf + sum(Ks.*(Ks+1)/2);
S.nnzAt    = nnz(Atf);
S.normb    = norm(full(bf));
S.c_nnz    = cnz;
S.feas     = (cnz==0);
% consistency: nvar must equal the cone dimension sossolve assumes
S.dim_ok   = (nvar == Kf + sum(double(Ks).^2));
end
