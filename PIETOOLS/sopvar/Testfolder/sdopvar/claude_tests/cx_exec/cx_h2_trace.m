function tr = cx_h2_trace(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TR = CX_H2_TRACE(X) is the trace of the finite-dimensional part of the
% square container X, as a 1x1 'dpvar': the sum of the diagonal entries of
% every diagonal block X.C{k,k} whose space k is R^n. It stands in for the
% executives' trace(Wm.P) (1-D), trace(tempObj.P) (1-D coercive) and
% trace(tempObj.R00) (2-D), so that the stock scalar lpi_ineq(prog,gam-TR)
% applies unchanged.
%
% Why not sdvar2dpvar (sopvar/Testfolder/converters): it indexes the
% constant row of matrix row r as (r-1)*nZd+1, but a dpvar holds nZd+1 rows
% per matrix row, so for an m x m block with m >= 2 the rows overlap. The
% trace needs only the diagonal, so it is summed here into a 1x1 dpvar
% directly and that converter is not used.
%
% Block layout (sdopvar.m header): an R^n -> R^n block has no variables, so
% its kernel is unvec(A + B'*d), m x m; entry (i,i) sits at vec position
% (i-1)*m+i in either vec order. A cdopvar holds ONE decision list X.Zd
% shared by every sdopvar block (cdopvar.m:39-40), which is the dpvar's
% dvarname list. With no R^n diagonal block the trace is the constant 0,
% returned as a dpvar so gam - TR stays a dpvar at fixed numeric gam (the
% stock 1-D non-coercive branch behaves so: trace of a 0x0 dpvar is a
% 1x1 dpvar).
%
% Cost: O(nnz of the diagonal columns of B); nothing densified along the
% decision axis.
%
% Initial coding MMP, 09/25/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isa(X,'cdopvar'),    Zd = X.Zd(:);   else,   Zd = cell(0,1);     end
q = numel(Zd);
a = 0;      b = sparse(q,1);
isR = ~any(X.space_out,2);                  % spaces with no variable: R^n
for k = find(isR(:)).'
    blk = X.C{k,k};
    if isempty(blk),    continue,   end
    if ~isempty(blk.vars.out) || ~isempty(blk.vars.in)
        error('cx_h2_trace:vars','Block (%d,%d) of an R^n space carries variables.',k,k)
    end
    m = blk.dims(1);
    idx = (0:m-1)*m + (1:m);                % diagonal positions in vec
    if isa(blk,'sdopvar')
        A = blk.params.A{1};    B = blk.params.B{1};
        if ~isequal(blk.Zd(:),Zd)
            error('cx_h2_trace:Zd','Block decision list differs from the container''s.')
        end
        if numel(A)==m*m,   a = a + full(sum(A(idx)));
        elseif any(A(:)),   error('cx_h2_trace:A','Nonzero scalar A shorthand in an m>1 block.')
        end                                 % [] or scalar 0: the zero shorthand
        if ~isempty(B), b = b + sum(B(:,idx),2);    end
    else                                    % fixed sopvar block: constants only
        C = blk.params{1};
        if ~isempty(C),  a = a + full(trace(C));   end
    end
end
tr = dpvar(sparse([a; b]),zeros(1,0),{},Zd,[1 1]);
end
