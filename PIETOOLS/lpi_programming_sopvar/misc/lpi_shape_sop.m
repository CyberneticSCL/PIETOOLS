function S = lpi_shape_sop(prog)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = LPI_SHAPE_SOP(PROG) the SDP shape of a (solved) LPI program as
% 'sossolve' hands it to the solver: S.ndv decision variables, S.Kf free
% block, S.Ks PSD block sizes (descending), S.m equality rows, S.nx primal
% length. Call it on the solved program: the 'ineq' slack blocks are added
% during the solve.
%
% Initial coding MMP, 10/08/2026 (the count of the test helper cx_shape,
% in the library so that the executives depend on the library only).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = struct('ndv',numel(prog.decvartable),'Kf',NaN,'Ks',[],'m',0,'nx',0);
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
S.nx = Kf + sum(S.Ks.^2);
end
