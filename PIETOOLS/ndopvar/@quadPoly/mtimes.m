function H = mtimes(F, G)
% H = MTIMES(F,G) computes the matrix product H of matrix-valued 'quadPoly'
% objects F and G
%
% INPUTS
% - F:  m x n 'quadPoly' or numeric object;
% - G:  n2 x p 'quadPoly' or numeric object. Unless F or G is scalar
%       (m=n=1 or n2=p=1), we must have n=n2;
%
% OUTPUTS
% - H:  m x p 'quadPoly' object representing the matrix product F*G; 
%
% NOTES
% - Left variables stay left (ns/Zs), right variables stay right (nt/Zt).
% As such, left variables in F cannot appear on right in G, or vice versa;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 06/25/2026: Initial coding;

% ---------- Deal with numeric * quadPoly case ----------
if isnumeric(F)
    H = leftNumericTimes(F, G);
    return;
elseif isa(F,'polynomial')
    if all(ismember(F.varname',[G.ns,G.nt]))
        F = polynomial2quadPoly(F,G.ns,G.nt);
    else
        error("Multiplication is not supported: 'polynomial' to 'quadPoly' conversion is ambiguous.")
    end
elseif ~isa(F,'quadPoly')
    error("Object must be either numeric or 'quadPoly' for multiplication.")
end
% ---------- Deal with quadPoly * numeric case ----------
if isnumeric(G)
    H = rightNumericTimes(F, G);
    return;
elseif isa(G,'polynomial')
    if all(ismember(G.varname',[F.ns,F.nt]))
        G = polynomial2quadPoly(G,F.ns,F.nt);
    else
        error("Multiplication is not supported: 'polynomial' to 'quadPoly' conversion is ambiguous.")
    end
elseif ~isa(G,'quadPoly')
    error("Object must be either numeric or 'quadPoly' for multiplication.")
end

% ---------- Deal with quadPoly * quadPoly case ----------
[m,n] = size(F);
[n2,p] = size(G);
if n ~= n2
    error('quadPoly:mtimes:dimMismatch', 'Inner dimensions must match.');
end

% Separately consider the case where F has only left monomials and G only
% right monomials
if isempty(F.nt) && isempty(G.ns)
    H = quadPoly(F.C*G.C, F.Zs, G.Zt, [m p], F.ns, G.nt);
    return
end

% Express objects in terms of same variables
if any(ismember(F.ns,G.nt)) || any(ismember(F.nt,G.ns))
    error("Multiplication is not supported: left-variables in one 'quadPoly' appear as right-variables in the other.")
end
[F,G] = common_vars(F,G);

% Take the product of the basis
[ZsP, sMapMats, nZs_arr] = productBasis1D(F.Zs, G.Zs);  % left variables: add exponents
[ZtP, tMapMats, nZt_arr] = productBasis1D(F.Zt, G.Zt);  % right variables
% Keep track of comulative dimension of monomial basis (for Kronecker
% product)
nnZs_arr = [cumprod(nZs_arr,'reverse'),1];  nZs = nnZs_arr(1);
nnZt_arr = [cumprod(nZt_arr,'reverse'),1];  nZt = nnZt_arr(1);

% Determine row index, column index, and left and right monomial indices 
% associated with all nonzero coefficients in F and G
nZs_F = cellfun(@(a) numel(a),F.Zs);
nZt_F = cellfun(@(a) numel(a),F.Zt);
[rsubs_F,csubs_F,vals_F] = find_idcs(F.C,[m,nZs_F],[n,nZt_F]);
nZs_G = cellfun(@(a) numel(a),G.Zs);
nZt_G = cellfun(@(a) numel(a),G.Zt);
[rsubs_G,csubs_G,vals_G] = find_idcs(G.C,[n,nZs_G],[p,nZt_G]);

% Take the product of all terms appearing in column i of F and row i of G
cidcs_F = csubs_F(:,1);
ridcs_G = rsubs_G(:,1);
rridcs_H = zeros(0,1);
ccidcs_H = zeros(0,1);
vals_H = zeros(0,1);
for i=1:n
    % Determine which of the nonzero coefficients appear in this row/column
    rtn_idcsF = cidcs_F==i;         nrtn_F = nnz(rtn_idcsF);
    rtn_idcsG = ridcs_G==i;         nrtn_G = nnz(rtn_idcsG);
    if nrtn_F==0 || nrtn_G==0
        continue
    end
    % Take all products of the coefficients in this row/column
    vals_i = vals_F(rtn_idcsF).*vals_G(rtn_idcsG)';
    vals_H = [vals_H; vals_i(:)];
    % Determine associated row and column indices in H
    ridcs_H_i = repmat(rsubs_F(rtn_idcsF,1),nrtn_G,1);
    cidcs_H_i = repelem(csubs_G(rtn_idcsG,1),nrtn_F,1);
    % Determine associated left and right monomial indices in F and G
    rsubs_F_i = repmat(rsubs_F(rtn_idcsF,2:end),nrtn_G,1);
    csubs_F_i = repmat(csubs_F(rtn_idcsF,2:end),nrtn_G,1);
    rsubs_G_i = repelem(rsubs_G(rtn_idcsG,2:end),nrtn_F,1);
    csubs_G_i = repelem(csubs_G(rtn_idcsG,2:end),nrtn_F,1);
    % Determine associated left and right monomial indices in H
    rsubs_H = pairToProdIndex(rsubs_F_i, rsubs_G_i, sMapMats);
    csubs_H = pairToProdIndex(csubs_F_i, csubs_G_i, tMapMats);
    % Convert to row and column indices in coefficient matrix for H
    rridcs_H_i = ([ridcs_H_i,rsubs_H]-1)*nnZs_arr' + 1;
    rridcs_H = [rridcs_H; rridcs_H_i];
    ccidcs_H_i = ([cidcs_H_i,csubs_H]-1)*nnZt_arr' + 1;
    ccidcs_H = [ccidcs_H; ccidcs_H_i];
end

% Declare the product H = F*G
C_H = sparse(rridcs_H,ccidcs_H,vals_H,m*nZs,p*nZt);
H = quadPoly(C_H, ZsP, ZtP, [m p], F.ns, F.nt);

end



%%
function H = leftNumericTimes(A, F)
% Perform multiplication of numeric A with quadPoly F
[q, m] = size(A);
if q==1 && m==1
    % Multiplication with scalar
    H = quadPoly(A*F.C, F.Zs, F.Zt, F.dim, F.ns, F.nt);
    return
elseif m ~= F.dim(1)
    error('quadPoly:mtimes:leftNumericDimMismatch', 'Left numeric matrix has wrong size.');
end

% Left multiplication: (A*F)(s,t) = A * F(s,t)
% Coeff update: Cnew = (kron(A, I_nZ)) * C
nZs = size(F.C,1)/m;
L = kron(sparse(A), speye(nZs));
Cnew = L * sparse(F.C);
H = quadPoly(Cnew, F.Zs, F.Zt, [q, F.dim(2)], F.ns, F.nt);

end


%%
function H = rightNumericTimes(F, B)
% Perform multiplication of quadPoly F with numeric B
[n, p] = size(B);
if n==1 && p==1
    % Multiplication with scalar
    H = quadPoly(F.C*B, F.Zs, F.Zt, F.dim, F.ns, F.nt);
    return
elseif n ~= F.dim(2)
    error('quadPoly:mtimes:rightNumericDimMismatch', 'Right numeric matrix has wrong size.');
end

% Right multiplication: (F*B)(s,t) = F(s,t) * B
% Coeff update: Cnew = C * kron(B, I_nZ)
nZt = size(F.C,2)/n;
R = kron(sparse(B), speye(nZt));
Cnew = sparse(F.C) * R;
H = quadPoly(Cnew, F.Zs, F.Zt, [F.dim(1), p], F.ns, F.nt);

end



%%
function [ZP, mapMats, nZP] = productBasis1D(ZA, ZB)
% Take the product of the monomial bases defined by ZA and ZB, by summing
% all degrees in ZA with all degrees in ZB
%
% Inputs:
%   ZA, ZB : 1×k cells, sorted-unique exponent vectors per variable
% Outputs:
%   ZP      : 1×k cells, sorted-unique sums per variable
%   mapMats : 1×k cells, mapMats{i} is (nAi×nBi) matrix giving index in ZP{i}
%   nZP     : 1 x k array of integers specifying the new number of exponents

k = numel(ZA);
ZP = cell(1,k);
mapMats = cell(1,k);
nZP = zeros(1,k);
for i = 1:k
    a = ZA{i}(:);
    b = ZB{i}(:);

    % All pairwise sums in this variable
    S = a + b.';                         % nAi × nBi
    z = unique(S(:), 'sorted');          % sorted-unique sums
    ZP{i} = z;
    nZP(i) = numel(z);

    % Map each sum back to its index in z
    [~, pos] = ismember(S(:), z);
    mapMats{i} = reshape(pos, size(S,1), size(S,2));
end

end



%%
function [rsubs,csubs,vals] = find_idcs(C,nr,nc)
% Find indices in monomial bases associated with each nonzero coefficient
% in the matrix C
%
% INPUTS
% - C:  (sparse) array of dimensions prod(nr) x prod(nc), representing
%       coefficient matrix acting on monomials 
%           ZL = ZL{1} o ... o ZL{m}    and     ZR = ZR{1} o ... o ZR{n}
%       for o the Kronekcer product
% - nr: 1 x m array specifying the number of elements in each monomial
%       basis ZL{i};
% - nc: 1 x n array specifying the number of elements in each monomial
%       basis ZR{i};
%
% OUTPUTS
% - rsubs:  k x m array specifying the monomial index in ZL associated with
%           each nonzero element of C;
% - csubs:  k x m array specifying the monomial index in ZL associated with
%           each nonzero element of C;
% - vals:   k x 1 array specifying the nonzero values appearing in C
%

% Find nonzero elements of C
[rridcs,ccidcs,vals] = find(C);
vals = vals(:);

% Split row indices into indices for different monomials
nnr = [cumprod(nr(2:end),'reverse'),1];     M = numel(nr);
rsubs = zeros(numel(rridcs),M);
rsubs(:,1) = ceil(rridcs(:)/nnr(1));
for i=1:M-1
    % Remove contribution of monomials in variable i
    rridcs = rridcs(:) - (rsubs(:,i)-1)*nnr(i);
    % Set indices associated with variable i+1
    rsubs(:,i+1) = ceil(rridcs/nnr(i+1));
end
% Split column indices into indices for different monomials
nnc = [cumprod(nc(2:end),'reverse'),1];     N = numel(nc);
csubs = zeros(numel(ccidcs),N);
csubs(:,1) = ceil(ccidcs(:)/nnc(1));
for i=1:N-1
    % Remove contribution of monomials in variable i
    ccidcs = ccidcs(:) - (csubs(:,i)-1)*nnc(i);
    % Set indices associated with variable i+1
    csubs(:,i+1) = ceil(ccidcs/nnc(i+1));
end

end



%%
function subsC = pairToProdIndex(subsA, subsB, mapMats)
% SUBSC = PAIRTOPRODINDEX(SUBSA,SUBSB,MAPMATS) takes the monomial indices 
% SUBSA and SUBSB of A and B, respectively, and determines associated 
% monomial indices in C, based on MAPMATS;
%
% INPTUS
% - subsA:  k x N array of integers, with element (i,j) specifying which 
%           monomial in variable j in A appears in term i;
% - subsB:  k x N array of integers, with element (i,j) specifying which 
%           monomial in variable j in B appears in term i;
% - mapMats:    1 x N cell with each element a dA_j x dB_j matrix of
%               integers, with element (p,q) indicating which monomial in
%               variable j in C corresponds with monomial p in A and
%               monomial q in B;
%
% OUTPUTS
% - subsC:  k x N array of integers, with element (i,j) specifying which 
%           monomial in variable j in C appears in term i;
%

subsC = zeros(size(subsA));
for i=1:size(subsC,2)
    mapMat_i = mapMats{i};
    lidcs = sub2ind(size(mapMat_i),subsA(:,i),subsB(:,i));
    subsC(:,i) = mapMat_i(lidcs);
end

end