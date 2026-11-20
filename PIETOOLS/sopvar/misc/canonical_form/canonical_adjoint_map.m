function [Tcell,ZLnew,ZRnew] = canonical_adjoint_map(vars,ZL,ZR,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [TCELL,ZLNEW,ZRNEW] = CANONICAL_ADJOINT_MAP(VARS,ZL,ZR,DIMS) returns the
% index maps taking each parameter of a 'sopvar' or 'sdopvar' operator to the
% corresponding parameter of its adjoint, keeping the multiplier parameters
% in canonical form.
%
% For an operator
%
%   (P x)(s) = sum_gam int_{s'} K_gam(s,s') I_gam(s-s') x(s') ds'
%   K_gam(s,s') = (I_m o ZL(s))^T unvec(C_gam) (I_n o ZR(s'))
%
% the adjoint is obtained by exchanging lower and upper integrals, 1->1,
% 2->3, 3->2, and transposing each kernel. In a direction carrying an
% integral that exchanges the roles of the primary and dummy variables, so
% the left and right monomial bases swap. In a direction carrying a
% multiplier it does NOT: delta(s_k-s_k') collapses the kernel onto the
% diagonal, so the parameter is a function of s_k alone,
%
%   M(s) = (I_m o ZL(s))^T unvec(C) (I_n o ZR(s)),
%
% and the adjoint of a multiplier is M(s)^T, the same polynomial in s with
% only the matrix indices transposed. The monomial degree therefore stays on
% the ZL side, and the ZR side stays at degree zero.
%
% This is what keeps the canonical form
%
%   gam_k == 1  ==>  the parameter has no dummy variable degree in
%                    direction k, that is, no content outside the columns
%                    whose ZR monomial has degree 0 in direction k
%
% invariant under taking adjoints. Applying only the plain vectorized
% transpose, as a naive implementation does, leaves the degrees on the ZR
% side and breaks it.
%
% INPUT
% - vars:   'struct' with fields 'in' and 'out', each a cellstr of variable
%           names, as stored by 'sopvar'/'sdopvar';
% - ZL:     cell array of column vectors of exponents, one per variable in
%           vars.out;
% - ZR:     cell array of column vectors of exponents, one per variable in
%           vars.in;
% - dims:   1x2 array [m,n] of row and column dimensions;
%
% OUTPUT
% - Tcell:  cell array of the same size as the parameter array. Tcell{k} is
%           the sparse 0/1 matrix with
%
%               vec(Cadj_{adj(k)}) = Tcell{k} * vec(C_k),
%
%           where adj(k) is the multi-index of k with 2 and 3 exchanged. Each
%           column holds at most one nonzero: entries that the canonical form
%           requires to be zero are dropped, so the input must be in
%           canonical form;
% - ZLnew:  left monomial basis of the adjoint, indexed by vars.in;
% - ZRnew:  right monomial basis of the adjoint, indexed by vars.out.
%
% NOTES
% Because ZL and ZR are properties of the object while the multi-index is a
% property of each parameter, a common direction generally carries a
% multiplier in some parameters and an integral in others. The adjoint bases
% must therefore serve both, and are formed as unions:
%
%   ZLnew{p} = union(ZR{p}, ZL{posL})    p a common variable of vars.in
%   ZLnew{p} = ZR{p}                     p an input-only variable
%   ZRnew{r} = union(ZL{r}, 0)           r a common variable of vars.out
%   ZRnew{r} = ZL{r}                     r an output-only variable
%
% For the usual case of one complete basis 0:deg shared by both sides, both
% unions are the identity and the bases are unchanged.
%
% See also IS_CANONICAL_MULTIPLIER, SOPVAR, SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - canonical_adjoint_map
%
% Copyright (C) 2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% MMP, 08/29/2026: Initial coding

% Variables, dimensions, and one basis per variable on each side
vin  = reshape(vars.in,1,[]);
vout = reshape(vars.out,1,[]);
Ni = numel(vin);        No = numel(vout);
m  = dims(1);           n  = dims(2);

ZL = reshape(ZL,1,[]);      ZR = reshape(ZR,1,[]);
if numel(ZL)~=No
    error("Left monomial basis should hold one entry per variable in 'vars.out'.")
end
if numel(ZR)~=Ni
    error("Right monomial basis should hold one entry per variable in 'vars.in'.")
end

% The parameter array is indexed by the variables common to the input and
% output spaces, in sorted order; those are the directions that can carry a
% multiplier.
vS3 = intersect(vin,vout);
n3  = numel(vS3);
[~,posL] = ismember(vS3,vout);      % where a common variable sits in ZL
[~,posR] = ismember(vS3,vin);       % where a common variable sits in ZR

nLold = cellfun(@numel,ZL);     NLold = prod([nLold,1]);
nRold = cellfun(@numel,ZR);     NRold = prod([nRold,1]);

% % % Bases of the adjoint. A common direction has to serve both the
% % % multiplier parameters, which keep their degree on the left, and the
% % % integral parameters, which swap sides.
% 'union' returns a ROW when both of its arguments are rows, and a 1x1 array
% counts as a row, so every result is forced back to a column: monomial
% bases are column vectors.
% Left basis of the adjoint, over vars.in
ZLnew = cell(1,Ni);
for p = 1:Ni
    t = find(posR==p,1);
    if isempty(t)
        ZLnew{p} = ZR{p}(:);
    else
        ZLnew{p} = union(ZR{p}(:),ZL{posL(t)}(:));
    end
    ZLnew{p} = ZLnew{p}(:);
end
% Right basis of the adjoint, over vars.out
ZRnew = cell(1,No);
for r = 1:No
    t = find(posL==r,1);
    if isempty(t)
        ZRnew{r} = ZL{r}(:);
    else
        ZRnew{r} = union(ZL{r}(:),0);
    end
    ZRnew{r} = ZRnew{r}(:);
end

nLnew = cellfun(@numel,ZLnew);      NLnew = prod([nLnew,1]);
nRnew = cellfun(@numel,ZRnew);      NRnew = prod([nRnew,1]);

% Monomials are ordered as ZL(x) = x_1^ZL{1} o ... o x_N^ZL{N}, so the first
% variable is the slowest index.
sLold = kron_strides(nLold);    sRold = kron_strides(nRold);
sLnew = kron_strides(nLnew);    sRnew = kron_strides(nRnew);

% Degree-to-index lookups, used to place a degree in the adjoint bases.
mapLnew = cell(1,Ni);
for p = 1:Ni
    mapLnew{p} = monomial_position_map(ZLnew{p});
end
mapRnew = cell(1,No);
for r = 1:No
    mapRnew{r} = monomial_position_map(ZRnew{r});
end

% Sizes of vec(C) before and after; one map per parameter cell
nrow_old = m*NLold;     ncol_old = n*NRold;     nC_old = nrow_old*ncol_old;
nrow_new = n*NLnew;     ncol_new = m*NRnew;     nC_new = nrow_new*ncol_new;

sz_C = [3*ones(1,n3),1];
Tcell = cell(sz_C);

for k = 1:3^n3
    gam = gamma_of_cell(k,n3);      % 1 x 0 when n3 = 0
    is_mult = (gam==1);

    % % % Only source columns whose ZR degree is zero in every multiplier
    % % % direction can carry content; the canonical form forbids the rest.
    okB = true(NRold,1);
    for t = find(is_mult)
        p = posR(t);
        v = (ZR{p}(:)==0);
        kk = true(1,1);
        for i = 1:Ni
            if i==p
                kk = kron(kk,v);
            else
                kk = kron(kk,true(nRold(i),1));
            end
        end
        okB = okB & kk;
    end
    bList = find(okB)-1;                % zero based, admissible ZR indices
    if isempty(bList)
        Tcell{k} = sparse(nC_new,nC_old);
        continue
    end

    % % % Enumerate the admissible source positions (i,aIdx,j,bIdx).
    [aG,iG,bG,jG] = ndgrid((0:NLold-1)',(0:m-1)',bList,(0:n-1)');
    aG = aG(:);     iG = iG(:);     bG = bG(:);     jG = jG(:);

    src = (jG*NRold + bG)*nrow_old + (iG*NLold + aG) + 1;

    alpha = kron_split(aG,sLold);       % zero based, per variable
    beta  = kron_split(bG,sRold);

    % % % Left index of the adjoint, over vars.in. A multiplier direction
    % % % keeps the degree it had on the left; every other direction takes
    % % % the degree it had on the right.
    aStar = zeros(numel(aG),1);
    for p = 1:Ni
        t = find(posR==p,1);
        if ~isempty(t) && is_mult(t)
            deg = ZL{posL(t)}(alpha(:,posL(t))+1);
        else
            deg = ZR{p}(beta(:,p)+1);
        end
        aStar = aStar + sLnew(p)*monomial_position(mapLnew{p},deg, ...
                                               "Adjoint monomial basis");
    end

    % % % Right index of the adjoint, over vars.out. A multiplier direction
    % % % is pinned at degree zero.
    bStar = zeros(numel(aG),1);
    for r = 1:No
        t = find(posL==r,1);
        if ~isempty(t) && is_mult(t)
            deg = zeros(numel(aG),1);
        else
            deg = ZL{r}(alpha(:,r)+1);
        end
        bStar = bStar + sRnew(r)*monomial_position(mapRnew{r},deg, ...
                                               "Adjoint monomial basis");
    end

    % Position in vec of the transposed coefficient matrix: row (j,aStar),
    % column (i,bStar)
    dest = (iG*NRnew + bStar)*nrow_new + (jG*NLnew + aStar) + 1;

    Tcell{k} = sparse(dest,src,1,nC_new,nC_old);
end

end
