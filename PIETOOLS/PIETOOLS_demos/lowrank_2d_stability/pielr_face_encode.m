function [W,res,info] = pielr_face_encode(Zop,Slist)                        % CC, 09/20/2026
% PIELR_FACE_ENCODE  Write known operators as columns of the Gram basis.
%
%   [W,res,info] = pielr_face_encode(Zop,Slist)
%
% WHY THIS EXISTS.  theory.pdf Section 2 reads the Gram face as the operator
% SQUARE ROOT written in the basis: with P = Z'*Q*Z and Q = Y*Y', each column
% of Y is one "square", i.e. the coefficient column of one ROW operator w'*Z.
% That reading is what predicts both measured halves of the rank law, and it
% is what would let V be WRITTEN DOWN from a known energy instead of searched
% (theory.pdf Section 4, "Physics-informed construction").  This routine is the
% read-out that makes the reading testable and, if it holds, usable: it maps an
% analytic operator to the Gram coordinates the face lives in.
%
% Zop is the monomial operator poslpivar_2d returns as its FOURTH output for
% the block concerned -- the opvar2d taking the state to the stacked basis
% column Z, so that the block's positive operator is exactly Zop'*Q*Zop.  Its
% rows are the Gram indices, in the SAME order the Gram vector uses.
%
% Slist is an opvar2d, or a cell of them, each a single-output ROW operator on
% the same state space.  For each, the w with w'*Zop = S is returned as a
% column of W.
%
% res(k) is the RELATIVE residual of that fit, measured on the stacked
% coefficient system.  res(k) ~ eps means S lies exactly in the basis image and
% W(:,k) is its unique coordinate column; res(k) = O(1) means the analytic
% operator is NOT representable at these degrees and W(:,k) is only a
% least-squares shadow -- which must be reported, never silently used, because
% a shadow column still produces a perfectly valid-looking face.
%
% The solution is unique whenever it exists: distinct basis rows carry distinct
% monomials in distinct R22 cells, so Zop has no left nullspace and the map
% w -> w'*Zop is injective.  (The map Q -> Zop'*Q*Zop that builds the OPERATOR
% may still have a nullspace; that is a separate question and this routine does
% not speak to it.)
%
% SCOPE.  R22 only: the 2-D-state block that every 2-D PDE stability LPI in
% this package produces (Top.dim = [0 0;0 0;0 0;n n]).  Anything outside R22 is
% an ERROR rather than a silent zero, since ignoring a populated cell would
% under-determine w and report a good residual for a wrong column.
%
% CC, 09/20/2026: initial version, for the candidate-M6 face read-out.

if ~iscell(Slist), Slist = {Slist}; end
K = numel(Slist);
N = size(Zop.R22{1,1},1);

% refuse anything with content outside R22 (see SCOPE)
chk = [{Zop},Slist(:)'];
for k = 1:numel(chk)
    pr = properties(chk{k});
    for p = 1:numel(pr)
        f = pr{p};
        if any(strcmp(f,{'I','var1','var2','dim','dimdependent','R22'})), continue, end
        v = chk{k}.(f);
        if ~iscell(v), v = {v}; end
        for c = 1:numel(v)
            assert(mxc_l(v{c})==0, ...
                'pielr_face_encode: operator %d has nonzero %s; only R22 is handled',k,f);
        end
    end
end

% stack one coefficient-matching system per R22 cell
A = [];  B = [];
for i = 1:3
    for j = 1:3
        Sc = cell(1,K);
        for k = 1:K, Sc{k} = Slist{k}.R22{i,j}; end
        [Aij,Bij] = cellsys_l(Zop.R22{i,j},Sc,N);
        A = [A;Aij];  B = [B;Bij];                                          %#ok<AGROW>
    end
end

W = A\B;
res = zeros(1,K);
for k = 1:K
    nb = norm(B(:,k));
    res(k) = norm(A*W(:,k)-B(:,k))/max(nb,eps);
end
info.A = A;  info.B = B;  info.N = N;  info.nrows = size(A,1);
info.rankA = rank(full(A));
end

% =========================================================================
function [A,b] = cellsys_l(Zc,Sc,N)
% One R22 cell: rows = the monomials appearing in EITHER side, columns = Gram
% indices.  A(:,k) is basis element k's coefficient vector, b(:,q) is target q.
K = numel(Sc);
[vZ,dZ,cZ] = pol3_l(Zc,N);
vS = cell(1,K);  dS = cell(1,K);  cS = cell(1,K);
for k = 1:K, [vS{k},dS{k},cS{k}] = pol3_l(Sc{k},1); end

vall = vZ(:);
for k = 1:K, vall = union(vall,vS{k}(:)); end
vall = vall(:);          % union of two ROW cellstrs returns a ROW (CLAUDE.md)

DZ = expand_l(dZ,vZ,vall);
DS = cell(1,K);
for k = 1:K, DS{k} = expand_l(dS{k},vS{k},vall); end

Dall = [DZ; vertcat(DS{:})];
if isempty(Dall), A = zeros(0,N); b = zeros(0,K); return, end
[~,~,ic] = unique(Dall,'rows');
nu = max(ic);
icZ = ic(1:size(DZ,1));
[r,c,v] = find(cZ);
A = sparse(icZ(r),c,v,nu,N);
b = zeros(nu,K);
off = size(DZ,1);
for k = 1:K
    ick = ic(off+(1:size(DS{k},1)));  off = off+size(DS{k},1);
    [r2,~,v2] = find(cS{k});
    if ~isempty(r2), b(:,k) = accumarray(ick(r2),v2,[nu,1]); end
end
end

% =========================================================================
function [vn,dm,C] = pol3_l(p,ncol)
% varname / degmat / (nterms x ncol) coefficient triple of one parameter cell,
% with the empty and plain-double cases normalised to a single constant term.
if isempty(p), vn = {}; dm = zeros(0,0); C = sparse(0,ncol); return, end
if isa(p,'double')
    vn = {};  dm = zeros(1,0);  C = sparse(reshape(p,1,[]));
    if size(C,2)~=ncol, C = sparse(1,ncol); C(1,:) = reshape(p,1,[]); end
    return
end
p  = polynomial(p);
vn = p.varname(:);
dm = p.degmat;
C  = p.coefficient;
if isempty(vn), dm = zeros(size(C,1),0); end
end

% =========================================================================
function D = expand_l(d,vfrom,vto)
% degmat over vfrom, re-expressed over the superset vto (missing vars: degree 0)
D = zeros(size(d,1),numel(vto));
if isempty(vfrom) || isempty(d), return, end
[tf,loc] = ismember(vfrom(:),vto(:));
assert(all(tf),'pielr_face_encode: variable missing from the union list');
D(:,loc) = d;
end

% =========================================================================
function m = mxc_l(X)
% max|coefficient| of one parameter cell (opcheck_2d's own mxc, kept local)
if isempty(X), m = 0; return, end
if     isa(X,'double'),     C = X;
elseif isa(X,'polynomial'), C = X.coefficient;
else,  m = NaN; return
end
if isempty(C), m = 0; else, m = full(max(abs(C(:)))); end
if isempty(m), m = 0; end
end
