function [params,ZL,ZR,changed] = canonicalize_multiplier(params,vars,ZL,ZR,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PARAMS,ZL,ZR,CHANGED] = CANONICALIZE_MULTIPLIER(PARAMS,VARS,ZL,ZR,DIMS)
% rewrites the parameters of a 'sopvar' or 'sdopvar' operator in canonical
% multiplier form, without changing the operator they represent.
%
% A parameter whose multi-index GAM has GAM(k)==1 carries a factor
% delta(s_k-s_k'), which identifies s_k with s_k'. A coefficient sitting at
% left degree a and right degree b in direction k therefore contributes the
% monomial s_k^(a+b), exactly as a coefficient at left degree a+b and right
% degree 0 would. The stored coefficients are consequently not determined by
% the operator. This routine removes that freedom by moving all of the degree
% onto the left,
%
%   (left a, right b)  -->  (left a+b, right 0),
%
% summing coefficients that land on the same monomial. The left basis is
% enlarged where it does not already hold the required degrees, and the right
% basis is enlarged to hold degree 0 if it does not.
%
% Parameters carrying an integral in a direction are untouched apart from
% being re-indexed into the enlarged bases: there the dummy variable is a
% genuine second variable and its degrees are not free.
%
% INPUT
% - params: parameters of the operator, either a cell array of coefficient
%           matrices ('sopvar') or a struct with fields 'A' and 'B' holding
%           cell arrays ('sdopvar');
% - vars:   'struct' with fields 'in' and 'out', each a cellstr;
% - ZL:     cell array of exponent vectors, one per variable in vars.out;
% - ZR:     cell array of exponent vectors, one per variable in vars.in;
% - dims:   1x2 array [m,n];
%
% OUTPUT
% - params: the canonical parameters;
% - ZL,ZR:  the monomial bases, enlarged if that was necessary;
% - changed: 'logical', false if the input was already canonical, in which
%           case nothing was rewritten.
%
% NOTES
% The cost when the input is already canonical, which is the normal case, is
% one sparse column count per parameter and nothing else. In particular it
% does not scale with the number of decision variables of an 'sdopvar'.
%
% See also IS_CANONICAL_MULTIPLIER, CANONICAL_ADJOINT_MAP, SOPVAR, SDOPVAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - canonicalize_multiplier
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

changed = false;

% % % Establish the structure, and bail out on anything not fully specified;
% % % validating the declaration is not this routine's job.
if ~isstruct(vars) || ~isfield(vars,'in') || ~isfield(vars,'out')
    return
end
vin  = reshape(vars.in,1,[]);
vout = reshape(vars.out,1,[]);
vS3  = intersect(vin,vout);
n3   = numel(vS3);

% Normalize the bases to a row of column vectors, the layout the class
% documents, so that adding a left and a right degree cannot silently expand
% into a matrix. This happens before the early returns below, so that an
% operator with no common variable is normalized just like any other.
ZL = reshape(ZL,1,[]);      ZR = reshape(ZR,1,[]);
for i = 1:numel(ZL),    ZL{i} = ZL{i}(:);    end
for i = 1:numel(ZR),    ZR{i} = ZR{i}(:);    end

if n3==0
    return              % no common direction, so no multiplier parameter
end
if numel(ZL)~=numel(vout) || numel(ZR)~=numel(vin) || numel(dims)<2
    return
end
[~,posL] = ismember(vS3,vout);
[~,posR] = ismember(vS3,vin);

is_sdop = isstruct(params);
if is_sdop
    if ~isfield(params,'A') || ~isfield(params,'B')
        return
    end
    ncell = numel(params.A);
else
    ncell = numel(params);
end
if ncell~=3^n3
    return
end

nL = cellfun(@numel,ZL);        NL = prod([nL,1]);
nR = cellfun(@numel,ZR);        NR = prod([nR,1]);
m = dims(1);        n = dims(2);
nrow = m*NL;        nC = nrow*(n*NR);

sL = strides_of(nL);        sR = strides_of(nR);
sz_C = [3*ones(1,n3),1];

% % % Pass 1. Locate the parameters that are not canonical, and record the
% % % left degrees that folding them will require.
need = cell(1,numel(vout));
for r = 1:numel(vout)
    need{r} = zeros(0,1);
end
bad = false(1,ncell);
need_zero = false(1,n3);
for k = 1:ncell
    idcs = cell(1,n3);
    [idcs{:}] = ind2sub(sz_C,k);
    gam = cell2mat(idcs);
    mult_dirs = find(gam==1);
    if isempty(mult_dirs)
        continue
    end
    % % % Positions the canonical form forbids: those whose right monomial
    % % % has a nonzero degree in some multiplier direction. Testing only
    % % % these costs O(numel(C)) and, when the parameter is canonical,
    % % % touches none of its nonzeros, so it is independent of the number
    % % % of decision variables.
    okB = true(NR,1);
    for t = mult_dirs
        p = posR(t);
        v = (ZR{p}(:)==0);
        kk = true(1,1);
        for i = 1:numel(vin)
            if i==p
                kk = kron(kk,v);
            else
                kk = kron(kk,true(nR(i),1));
            end
        end
        okB = okB & kk;
    end
    bBad = find(~okB)-1;
    if isempty(bBad) || ~has_content(params,k,is_sdop)
        continue                % nothing forbidden, or nothing stored at all
    end
    [bG,jG,rG] = ndgrid(bBad,(0:n-1)',(0:nrow-1)');
    linBad = (jG(:)*NR + bG(:))*nrow + rG(:) + 1;

    lin = linBad(forbidden_content(params,k,is_sdop,nC,linBad));
    if isempty(lin)
        continue
    end
    bad(k) = true;
    need_zero(mult_dirs) = true;
    row = mod(lin-1,nrow);          col = floor((lin-1)/nrow);
    alpha = split_index(mod(row,NL),sL);
    beta  = split_index(mod(col,NR),sR);
    for t = mult_dirs
        r = posL(t);
        tot = ZL{r}(alpha(:,r)+1) + ZR{posR(t)}(beta(:,posR(t))+1);
        need{r} = [need{r}; tot(:)];
    end
end
if ~any(bad)
    return
end
changed = true;

% % % Pass 2. Enlarge the bases. The left basis has to hold the folded
% % % degrees as well as the ones the integral parameters already use; the
% % % right basis has to hold degree 0.
% 'union' returns a ROW when both of its arguments are rows, and a 1x1 array
% counts as a row, so every result is forced back to a column. A row basis
% is not merely untidy: 'UnionBasisMonomials' builds its degree matrix from
% the orientation, and 'plus', 'minus' and 'eq' then fail outright.
ZLnew = ZL;     ZRnew = ZR;
for r = 1:numel(vout)
    if ~isempty(need{r})
        ZLnew{r} = union(ZL{r}(:),unique(need{r}));
    else
        ZLnew{r} = ZL{r};
    end
    ZLnew{r} = ZLnew{r}(:);
end
for t = find(need_zero)
    ZRnew{posR(t)} = union(ZR{posR(t)}(:),0);
end
for p = 1:numel(vin)
    ZRnew{p} = ZRnew{p}(:);
end

nLn = cellfun(@numel,ZLnew);        NLn = prod([nLn,1]);
nRn = cellfun(@numel,ZRnew);        NRn = prod([nRn,1]);
sLn = strides_of(nLn);              sRn = strides_of(nRn);
nrow_n = m*NLn;     nC_n = nrow_n*(n*NRn);

mapL = cell(1,numel(vout));
for r = 1:numel(vout),  mapL{r} = degree_lookup(ZLnew{r});   end
mapR = cell(1,numel(vin));
for p = 1:numel(vin),   mapR{p} = degree_lookup(ZRnew{p});   end

% % % Pass 3. Rebuild every parameter in the new bases, folding the
% % % multiplier directions on the way. Every parameter has to be rebuilt,
% % % not just the offending ones, because enlarging a basis moves the
% % % positions of all the others.
%
% A parameter whose size disagrees with the declared bases cannot be
% rebuilt, and silently emptying it would discard content, so refuse. Pass 1
% has already found a well formed parameter that needs folding, so reaching
% this with an inconsistent sibling means the parameter array is malformed.
for k = 1:ncell
    if is_sdop
        ok = (isempty(params.A{k}) || numel(params.A{k})==nC) ...
             && (isempty(params.B{k}) || size(params.B{k},2)==nC);
    else
        ok = isempty(params{k}) || numel(params{k})==nC;
    end
    if ~ok
        error("Parameter "+num2str(k)+" does not have "+num2str(nC)+" entries, as "...
              +"the declared dimensions and monomial bases require, so it cannot be "...
              +"rewritten in canonical form.")
    end
end
for k = 1:ncell
    idcs = cell(1,n3);
    [idcs{:}] = ind2sub(sz_C,k);
    gam = cell2mat(idcs);
    is_mult = (gam==1);

    lin = nonzero_positions(params,k,is_sdop,nC);
    if isempty(lin)
        params = set_param(params,k,is_sdop,[],nC_n,nrow_n);
        continue
    end
    row = mod(lin-1,nrow);      col = floor((lin-1)/nrow);
    i0 = floor(row/NL);         aIdx = mod(row,NL);
    j0 = floor(col/NR);         bIdx = mod(col,NR);
    alpha = split_index(aIdx,sL);
    beta  = split_index(bIdx,sR);

    aStar = zeros(numel(lin),1);
    for r = 1:numel(vout)
        t = find(posL==r,1);
        deg = ZL{r}(alpha(:,r)+1);
        if ~isempty(t) && is_mult(t)
            deg = deg + ZR{posR(t)}(beta(:,posR(t))+1);
        end
        aStar = aStar + sLn(r)*lookup(mapL{r},deg);
    end
    bStar = zeros(numel(lin),1);
    for p = 1:numel(vin)
        t = find(posR==p,1);
        if ~isempty(t) && is_mult(t)
            deg = zeros(numel(lin),1);
        else
            deg = ZR{p}(beta(:,p)+1);
        end
        bStar = bStar + sRn(p)*lookup(mapR{p},deg);
    end

    dest = (j0*NRn + bStar)*nrow_n + (i0*NLn + aStar) + 1;
    T = sparse(dest,lin,1,nC_n,nC);
    params = set_param(params,k,is_sdop,T,nC_n,nrow_n);
end

ZL = ZLnew;     ZR = ZRnew;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = has_content(params,k,is_sdop)
% Whether parameter k stores anything at all. Tested before the position
% enumeration below, so that an empty or all-zero parameter costs O(1)
% instead of one index per coefficient.

if is_sdop
    tf = nnz(params.A{k})>0 || nnz(params.B{k})>0;
else
    tf = nnz(params{k})>0;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = forbidden_content(params,k,is_sdop,nC,linBad)
% Which of the given positions of parameter k carry any content. Only the
% listed positions are inspected, so a canonical parameter costs
% numel(linBad) and none of its own nonzeros.

% Indexing a vector returns the orientation of the vector it indexes, not of
% the index, so every result is reshaped to a column before combining.
tf = false(numel(linBad),1);
if is_sdop
    Ak = params.A{k};       Bk = params.B{k};
    if ~isempty(Ak) && numel(Ak)==nC
        tf = tf | full(reshape(Ak(linBad)~=0,[],1));
    end
    if ~isempty(Bk) && size(Bk,2)==nC
        tf = tf | full(reshape(any(Bk(:,linBad),1),[],1));
    end
else
    Ck = params{k};
    if ~isempty(Ck) && numel(Ck)==nC
        tf = tf | full(reshape(Ck(linBad)~=0,[],1));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lin = nonzero_positions(params,k,is_sdop,nC)
% Linear positions within vec(C_k) that carry any content, whether through
% the constant part or through a decision variable.

if is_sdop
    Ak = params.A{k};       Bk = params.B{k};
    la = [];    lb = [];
    if ~isempty(Ak) && numel(Ak)==nC
        la = find(Ak);
    end
    if ~isempty(Bk) && size(Bk,2)==nC
        lb = find(any(Bk,1));
    end
    lin = unique([la(:); lb(:)]);
else
    Ck = params{k};
    if isempty(Ck) || numel(Ck)~=nC
        lin = [];
    else
        lin = find(Ck(:));
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = set_param(params,k,is_sdop,T,nC_n,nrow_n)
% Apply the index map to parameter k. A parameter with no content is resized
% to the new bases, or left empty if that is how it was stored, since callers
% treat an empty parameter as zero.

if is_sdop
    if isempty(T)
        if isempty(params.A{k}),    params.A{k} = [];
        else,                       params.A{k} = sparse(nC_n,1);
        end
        if isempty(params.B{k}),    params.B{k} = [];
        else,                       params.B{k} = sparse(size(params.B{k},1),nC_n);
        end
    else
        Ak = params.A{k};
        if isempty(Ak),     params.A{k} = sparse(nC_n,1);
        else,               params.A{k} = T*sparse(Ak(:));
        end
        Bk = params.B{k};
        if isempty(Bk),     params.B{k} = [];
        else,               params.B{k} = Bk*T.';
        end
    end
else
    if isempty(T)
        if isempty(params{k}),      params{k} = [];
        else,                       params{k} = sparse(nrow_n,nC_n/nrow_n);
        end
    else
        params{k} = reshape(T*sparse(params{k}(:)),nrow_n,nC_n/nrow_n);
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = strides_of(nvec)
% Stride of each variable in kron(Z{1},...,Z{N}), first variable slowest.

N = numel(nvec);
s = ones(1,N);
for k = N-1:-1:1
    s(k) = s(k+1)*nvec(k+1);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = split_index(idx,stride)
% Split a zero-based monomial index into zero-based per-variable indices.

N = numel(stride);
a = zeros(numel(idx),N);
rem = idx(:);
for k = 1:N
    a(:,k) = floor(rem/stride(k));
    rem = rem - a(:,k)*stride(k);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function map = degree_lookup(Zvec)
% Position of a degree within a basis, as a zero-based index.

d = double(Zvec(:));
if isempty(d)
    map = zeros(0,1);
    return
end
map = -ones(max(d)+1,1);
map(d+1) = 0:numel(d)-1;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = lookup(map,deg)
% Position of each requested degree in a basis, as a zero-based index.

deg = double(deg(:));
out = -ones(numel(deg),1);
in_range = deg>=0 & deg==round(deg) & deg+1<=numel(map);
out(in_range) = map(deg(in_range)+1);
bad = find(out<0,1);
if ~isempty(bad)
    error("Monomial basis does not contain degree "+num2str(deg(bad))...
          +"; the bases were not enlarged correctly.")
end

end
