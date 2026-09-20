function P = pielr_evalgram(prog,P)
% PIELR_EVALGRAM -- substitute a solved decision vector into a dopvar2d,
% directly, with sparse algebra only.
%
%   P = pielr_evalgram(prog,P)
%
% Returns the opvar2d obtained by evaluating every dpvar parameter cell of
% the dopvar2d P at the decision values prog.solinfo.RRx (ordered as
% prog.decvartable) -- the same operator lpigetsol/getsol_lpivar_2d returns,
% without going through sosgetsol.
%
% WHY THIS EXISTS.  On PIETOOLS master (998a68cc) sosgetsol handles a dpvar
% by dpvar2poly + combine, which expands every decision variable into a
% polynomial variable: on this package's negativity operator (~277k dvars,
% heat n=1) getsol_lpivar_2d did not finish in 35 minutes there, vs ~1 s on
% the development branch.  opcheck_2d / opcheck_pos2d were the package's
% only getsol dependence, so verification substitutes here instead and the
% package runs on both trees.
%
% LAYOUT (verified against the live class 2026-08-28; the dpvar classdef
% header claims the opposite and is WRONG): dpvar(C,degmat,varname,dvarname,
% [m,n]) stores C with ROWS matrix-row outer, decision-variable inner -- one
% block of nd+1 consecutive rows per matrix row, constant term FIRST, i.e. a
% left factor kron(I_m,[1 d']) -- and COLUMNS matrix-column outer, monomial
% inner.  Evaluation at d is therefore ONE sparse product per cell,
%     Cnew = kron(speye(m),sparse([1;d]')) * C          (m x n*nZ)
% followed by an index remap into polynomial's T x (m*n) column-major
% coefficient layout.  Nothing is ever full()'d on a decision-variable
% dimension (CLAUDE.md performance rules; rows of C are the ~q axis).
%
% Each cell's dvarname (a subset of the program's dvars, in its own order)
% is mapped to the global vector BY NAME against prog.decvartable, whose
% positions index RRx directly -- the same mapping sosgetsol's own dpvar
% path uses.  A dvar missing from decvartable is an ERROR, never a symbolic
% leftover: in this package every dvar comes from the program, so a miss
% means the wrong program/operator pairing, and hiding it would corrupt the
% operator-level verification this feeds.
%
% INPUT
%   prog - LPI/SOS program struct; only .decvartable (cellstr) and
%          .solinfo.RRx (decision vector, decvartable order) are used
%   P    - dopvar2d to evaluate; an opvar2d passes through unchanged
% OUTPUT
%   P    - opvar2d with every cell evaluated (polynomial or double)
%
% The 36-leaf enumeration below is getsol_lpivar_2d's own; if the class ever
% grows a parameter field, opcheck's 16-field/36-leaf assert on this
% function's output catches it.
%
% CC, 09/19/2026: initial version.

if isa(P,'opvar2d')
    return                              % nothing to substitute
end
if ~isa(P,'dopvar2d')
    error('pielr_evalgram:input','input must be a dopvar2d (or opvar2d) object.');
end
dtab = prog.decvartable;
RRx  = prog.solinfo.RRx;

for f = {'R00','R0x','R0y','R02','Rx0','Rxy','Ry0','Ryx','R20'}
    P.(f{:}) = evaldp(P.(f{:}),dtab,RRx);
end
for i = 1:3
    P.Rxx{i,1} = evaldp(P.Rxx{i,1},dtab,RRx);
    P.Rx2{i,1} = evaldp(P.Rx2{i,1},dtab,RRx);
    P.R2x{i,1} = evaldp(P.R2x{i,1},dtab,RRx);
    P.Ryy{1,i} = evaldp(P.Ryy{1,i},dtab,RRx);
    P.Ry2{1,i} = evaldp(P.Ry2{1,i},dtab,RRx);
    P.R2y{1,i} = evaldp(P.R2y{1,i},dtab,RRx);
    for j = 1:3
        P.R22{i,j} = evaldp(P.R22{i,j},dtab,RRx);
    end
end
P = opvar2d(P);                         % dopvar2d with polynomial cells -> opvar2d
end

% =========================================================================
function p = evaldp(V,dtab,RRx)
% Evaluate ONE parameter cell.  Only dpvar cells carry decision variables;
% empty/double/polynomial cells pass through, as in getsol_lpivar_2d.
if isempty(V) || isa(V,'double') || isa(V,'polynomial')
    p = V;  return
end
m  = V.matdim(1);   n  = V.matdim(2);
dn = V.dvarname;    nd = numel(dn);
dm = V.degmat;      vn = V.varname;
nZ = size(dm,1);
C  = V.C;
assert(size(C,1)==m*(nd+1) && size(C,2)==n*nZ, ...
    'pielr_evalgram: cell C is %dx%d, expected %dx%d from [m n nd nZ]=[%d %d %d %d]', ...
    size(C,1),size(C,2),m*(nd+1),n*nZ,m,n,nd,nZ);
if nd > 0
    [tf,loc] = ismember(dn,dtab);
    assert(all(tf), ...
        'pielr_evalgram: %d decision variable(s) of this cell are not in prog.decvartable',sum(~tf));
    d = full(RRx(loc));  d = d(:);
else
    d = zeros(0,1);
end
% row i of the evaluated coefficient = (const row) + d'*(dvar rows) of block i
Cnew = kron(speye(m),sparse([1;d]'))*C;     % m x n*nZ, stays sparse throughout
% remap to polynomial layout: entry (t,(j-1)*nZ+... ) -> row t, column i+m*(j-1)
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
