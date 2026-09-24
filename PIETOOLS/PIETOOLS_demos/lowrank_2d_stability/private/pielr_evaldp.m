function p = pielr_evaldp(V,dtab,RRx)                                       % CC, 09/23/2026
% PIELR_EVALDP  Evaluate ONE dpvar parameter cell at a decision vector.
%
% Body moved VERBATIM from the local function evaldp inside pielr_evalgram.m
% (CC, 09/19/2026), which is now one of two real callers: the 1-D evaluator
% pielr_evalop needs the identical per-cell arithmetic over the 6 leaves of a
% dopvar, and duplicating it would let the two copies drift -- exactly the
% failure this package already hit with the forked tests_1d core.
%
% Only dpvar cells carry decision variables; empty/double/polynomial cells
% pass through, as in getsol_lpivar_2d.
%
% LAYOUT (verified against the live class 2026-08-28; the dpvar classdef
% header claims the opposite and is WRONG): dpvar stores C with ROWS
% matrix-row outer, decision-variable inner -- one block of nd+1 consecutive
% rows per matrix row, constant term FIRST, i.e. a left factor
% kron(I_m,[1 d']) -- and COLUMNS matrix-column outer, monomial inner.
% Nothing is ever full()'d on the decision-variable axis: rows of C are the
% ~q dimension (CLAUDE.md section 2).

if isempty(V) || isa(V,'double') || isa(V,'polynomial')
    p = V;  return
end
m  = V.matdim(1);   n  = V.matdim(2);
dn = V.dvarname;    nd = numel(dn);
dm = V.degmat;      vn = V.varname;
nZ = size(dm,1);
C  = V.C;
assert(size(C,1)==m*(nd+1) && size(C,2)==n*nZ, ...
    'pielr_evaldp: cell C is %dx%d, expected %dx%d from [m n nd nZ]=[%d %d %d %d]', ...
    size(C,1),size(C,2),m*(nd+1),n*nZ,m,n,nd,nZ);
if nd > 0
    [tf,loc] = ismember(dn,dtab);
    assert(all(tf), ...
        'pielr_evaldp: %d decision variable(s) of this cell are not in prog.decvartable',sum(~tf));
    d = full(RRx(loc));  d = d(:);
else
    d = zeros(0,1);
end
% row i of the evaluated coefficient = (const row) + d'*(dvar rows) of block i
Cnew = kron(speye(m),sparse([1;d]'))*C;     % m x n*nZ, stays sparse throughout
% remap to polynomial layout: entry (t,(j-1)*nZ+...) -> row t, column i+m*(j-1)
[ri,ci,vv] = find(Cnew);
jj = floor((ci-1)/nZ);                      % matrix column minus 1
tt = ci - jj*nZ;                            % monomial index
if nZ==0 || isempty(vn)
    % no polynomial part: plain matrix (sparse() sums any duplicate indices)
    p = reshape(full(sparse(ri+m*jj,ones(size(ri)),vv,m*n,1)),m,n);
else
    p = combine(polynomial(sparse(tt,ri+m*jj,vv,nZ,m*n),dm,vn,[m,n]));
end
end
