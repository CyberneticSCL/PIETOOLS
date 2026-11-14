function S = cx_shape(prog)
% S = CX_SHAPE(PROG) is the SDP a (solved) PIETOOLS program handed the
% solver: S.ndv decision variables (numel(prog.decvartable)), S.Kf free
% block, S.Ks PSD block sizes (sorted descending, for multiset comparison),
% S.m equality rows, S.nvar primal length. Reproduces sossolve's own count
% (processvars is local to sossolve.m, so cannot be called), as
% PIETOOLS_demos/cuadmm/private/sdpshape.m does; call it on the SOLVED
% program, since sossolve adds the 'ineq' slack blocks during the solve.
%
% Initial coding MMP, 09/25/2026
S = struct('ndv',numel(prog.decvartable),'Kf',NaN,'Ks',[],'m',0,'nvar',0);
Kf = prog.var.idx{1}-1;
Ks = [];
for i = 1:prog.var.num
    sz = prog.var.idx{i+1}-prog.var.idx{i};
    switch prog.var.type{i}
        case 'poly',    Kf = Kf + sz;
        case 'sos',     Ks(end+1) = sqrt(sz);                               %#ok<AGROW>
    end
end
if isfield(prog,'extravar') && isfield(prog.extravar,'num')
    for i = 1:prog.extravar.num
        Ks(end+1) = sqrt(prog.extravar.idx{i+1}-prog.extravar.idx{i});    %#ok<AGROW>
    end
end
m = 0;
for i = 1:prog.expr.num,    m = m + size(prog.expr.At{i},2);    end
S.Kf = Kf;  S.Ks = sort(round(Ks),'descend');   S.m = m;
S.nvar = Kf + sum(S.Ks.^2);
end
