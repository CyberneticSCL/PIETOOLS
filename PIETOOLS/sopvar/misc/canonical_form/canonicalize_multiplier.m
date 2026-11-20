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
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

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

% % % Establish the structure; a declaration that is not fully specified
% % % is returned unchanged.
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
% Positions of the shared variables in ZL and in ZR
[~,posL] = ismember(vS3,vout);
[~,posR] = ismember(vS3,vin);

% 'sdopvar' parameters are a struct of A and B cells; one cell per
% multi-index
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

% Basis sizes and the size of vec(C)
nL = cellfun(@numel,ZL);        NL = prod([nL,1]);
nR = cellfun(@numel,ZR);        NR = prod([nR,1]);
m = dims(1);        n = dims(2);
nrow = m*NL;        nC = nrow*(n*NR);

% Strides of the kron(Z{1},...,Z{N}) bases, first variable slowest
sL = kron_strides(nL);      sR = kron_strides(nR);

% Pass 1 detects; passes 2 and 3 (enlarge the bases, rebuild the
% parameters) run only when it finds a non-canonical parameter, which none
% of the 3788 constructions of the 1-D, 2-D and 3-D executive program
% builds met (measured). Pass 1 runs on every construction, so it takes the
% sizes it reads as arguments (a struct built and unpacked costs 7.3 us per
% call, arguments 0.03 us); the rewrite passes take theirs as the struct
% 'sz', built on that path only.
[bad,need,need_zero] = find_noncanonical(params,ZL,ZR,vin,vout,n3,posL, ...
    posR,is_sdop,ncell,nR,NL,NR,n,nrow,nC,sL,sR);
if ~any(bad)
    return
end
changed = true;
sz = struct('vin',{vin},'vout',{vout},'n3',n3,'posL',posL,'posR',posR, ...
            'is_sdop',is_sdop,'ncell',ncell,'NL',NL,'NR',NR,'m',m, ...
            'n',n,'nrow',nrow,'nC',nC,'sL',sL,'sR',sR);
[ZLnew,ZRnew] = enlarge_bases(ZL,ZR,need,need_zero,sz);
params = rebuild_params(params,ZL,ZR,ZLnew,ZRnew,sz);
ZL = ZLnew;     ZR = ZRnew;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bad,need,need_zero] = find_noncanonical(params,ZL,ZR,vin,vout, ...
    n3,posL,posR,is_sdop,ncell,nR,NL,NR,n,nrow,nC,sL,sR)
% % % Pass 1. Locate the parameters that are not canonical, and record the
% % % left degrees that folding them will require.
% INPUTS  params, ZL, ZR, vin, vout, and the sizes of the primary function,
%         under the same names.
% OUTPUT  bad(k): parameter k has content the canonical form forbids;
%         need{r}: the left degrees of output variable r that folding
%         needs; need_zero(t): ZR needs degree 0 in shared direction t.
% Reads params, writes nothing.
need = cell(1,numel(vout));
for r = 1:numel(vout)
    need{r} = zeros(0,1);
end
bad = false(1,ncell);
need_zero = false(1,n3);
% maskR{t}: right monomials of degree 0 in shared direction t. Depends on
% (ZR,t) only, so built on first use and reused by every later cell.
maskR = cell(1,n3);
for k = 1:ncell
    % Empty cell first: nothing stored, so nothing to fold. It then costs
    % one nnz per stored array, not the multi-index, the O(NR) mask and one
    % index per coefficient.
    if is_sdop
        tf = stores(params.A{k}) || stores(params.B{k});
    else
        tf = stores(params{k});
    end
    if ~tf
        continue
    end
    gam = gamma_of_cell(k,n3);      % multi-index of cell k, labels 1/2/3
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
        if ~isempty(maskR{t})
            okB = okB & maskR{t};
            continue
        end
        % Mask of direction t on first use: degree 0 in variable p, any
        % degree in the others, over the kron basis
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
        maskR{t} = kk;
        okB = okB & kk;
    end
    bBad = find(~okB)-1;
    if isempty(bBad)                                % content tested above
        continue                % nothing forbidden
    end
    % Linear positions in vec(C) of the forbidden columns, and those that
    % hold content
    [bG,jG,rG] = ndgrid(bBad,(0:n-1)',(0:nrow-1)');
    linBad = (jG(:)*NR + bG(:))*nrow + rG(:) + 1;

    lin = linBad(forbidden_content(params,k,is_sdop,nC,linBad));
    if isempty(lin)
        continue
    end
    % Record the folded left degrees a+b of every multiplier direction
    bad(k) = true;
    need_zero(mult_dirs) = true;
    row = mod(lin-1,nrow);          col = floor((lin-1)/nrow);
    alpha = kron_split(mod(row,NL),sL);
    beta  = kron_split(mod(col,NR),sR);
    for t = mult_dirs
        r = posL(t);
        tot = ZL{r}(alpha(:,r)+1) + ZR{posR(t)}(beta(:,posR(t))+1);
        need{r} = [need{r}; tot(:)];
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ZLnew,ZRnew] = enlarge_bases(ZL,ZR,need,need_zero,sz)
% % % Pass 2. Enlarge the bases. The left basis has to hold the folded
% % % degrees as well as the ones the integral parameters already use; the
% % % right basis has to hold degree 0.
% INPUTS  ZL, ZR: the bases; need, need_zero: as returned by pass 1; sz:
%         'struct' of the sizes of the primary function.
% OUTPUT  ZLnew, ZRnew: the enlarged bases.
% 'union' returns a ROW when both of its arguments are rows, and a 1x1 array
% counts as a row, so every result is forced back to a column: monomial
% bases are column vectors.
vin = sz.vin;   vout = sz.vout;     posR = sz.posR;
ZLnew = ZL;     ZRnew = ZR;
for r = 1:numel(vout)
    if ~isempty(need{r})
        ZLnew{r} = union(ZL{r}(:),unique(need{r}));
    else
        ZLnew{r} = ZL{r};
    end
    ZLnew{r} = ZLnew{r}(:);
end
% Degree 0 in the right basis of each shared direction that needs it
for t = find(need_zero)
    ZRnew{posR(t)} = union(ZR{posR(t)}(:),0);
end
for p = 1:numel(vin)
    ZRnew{p} = ZRnew{p}(:);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = rebuild_params(params,ZL,ZR,ZLnew,ZRnew,sz)
% % % Pass 3. Rebuild every parameter in the new bases, folding the
% % % multiplier directions on the way. Every parameter has to be rebuilt,
% % % not just the offending ones, because enlarging a basis moves the
% % % positions of all the others.
% INPUTS  params: the parameters; ZL, ZR: the old bases; ZLnew, ZRnew: the
%         enlarged bases; sz: 'struct' of the sizes of the primary function.
% OUTPUT  params: the canonical parameters on the new bases.
%
% A parameter whose size disagrees with the declared bases cannot be
% rebuilt, and silently emptying it would discard content, so refuse. Pass 1
% has already found a well formed parameter that needs folding, so reaching
% this with an inconsistent sibling means the parameter array is malformed.
% Called as params = rebuild_params(params,...), the form MATLAB can update
% in place; each cell is replaced whole.
vin = sz.vin;   vout = sz.vout;     n3 = sz.n3;     posL = sz.posL;
posR = sz.posR;     is_sdop = sz.is_sdop;   ncell = sz.ncell;   m = sz.m;
n = sz.n;   NL = sz.NL;     NR = sz.NR;     nrow = sz.nrow;     nC = sz.nC;
sL = sz.sL;     sR = sz.sR;
% Sizes, strides and degree lookups of the new bases
nLn = cellfun(@numel,ZLnew);        NLn = prod([nLn,1]);
nRn = cellfun(@numel,ZRnew);        NRn = prod([nRn,1]);
sLn = kron_strides(nLn);            sRn = kron_strides(nRn);
nrow_n = m*NLn;     nC_n = nrow_n*(n*NRn);

mapL = cell(1,numel(vout));
for r = 1:numel(vout),  mapL{r} = monomial_position_map(ZLnew{r});   end
mapR = cell(1,numel(vin));
for p = 1:numel(vin),   mapR{p} = monomial_position_map(ZRnew{p});   end

% Every parameter must have the declared size
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
% Rebuild each parameter from its nonzero positions
for k = 1:ncell
    gam = gamma_of_cell(k,n3);
    is_mult = (gam==1);

    lin = nonzero_positions(params,k,is_sdop,nC);
    if isempty(lin)
        params = set_param(params,k,is_sdop,[],nC_n,nrow_n);
        continue
    end
    % Matrix entry (i0,j0) and per-variable monomial positions of each
    row = mod(lin-1,nrow);      col = floor((lin-1)/nrow);
    i0 = floor(row/NL);         aIdx = mod(row,NL);
    j0 = floor(col/NR);         bIdx = mod(col,NR);
    alpha = kron_split(aIdx,sL);
    beta  = kron_split(bIdx,sR);

    % New left index: a multiplier direction adds the right degree to the
    % left one
    aStar = zeros(numel(lin),1);
    for r = 1:numel(vout)
        t = find(posL==r,1);
        deg = ZL{r}(alpha(:,r)+1);
        if ~isempty(t) && is_mult(t)
            deg = deg + ZR{posR(t)}(beta(:,posR(t))+1);
        end
        aStar = aStar + sLn(r)*monomial_position(mapL{r},deg);
    end
    % New right index: degree zero in a multiplier direction
    bStar = zeros(numel(lin),1);
    for p = 1:numel(vin)
        t = find(posR==p,1);
        if ~isempty(t) && is_mult(t)
            deg = zeros(numel(lin),1);
        else
            deg = ZR{p}(beta(:,p)+1);
        end
        bStar = bStar + sRn(p)*monomial_position(mapR{p},deg);
    end

    % Index map into vec(C) on the new bases; folded entries sum
    dest = (j0*NRn + bStar)*nrow_n + (i0*NLn + aStar) + 1;
    T = sparse(dest,lin,1,nC_n,nC);
    params = set_param(params,k,is_sdop,T,nC_n,nrow_n);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = stores(X)
% TF = STORES(X) is true when X has a nonzero entry. nnz is O(1) for sparse
% storage but scans a full array; find(X,1) stops at the first nonzero.
% Pass 1 calls this for every cell, so a full q x nC B would otherwise cost
% O(q*nC) per cell.
% INPUT   X: numeric array, sparse or full.
% OUTPUT  tf: logical.
if issparse(X),     tf = nnz(X)>0;
else,               tf = ~isempty(find(X,1));
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = forbidden_content(params,k,is_sdop,nC,linBad)
% TF = FORBIDDEN_CONTENT(PARAMS,K,IS_SDOP,NC,LINBAD) tells which of the
% given positions of parameter k carry any content. Only the listed
% positions are inspected, so a canonical parameter costs numel(linBad) and
% none of its own nonzeros.
% INPUTS  params: the parameters; k: cell index; is_sdop: true for 'sdopvar'
%         parameters (A and B); nC: numel of vec(C); linBad: positions in
%         vec(C).
% OUTPUT  tf: numel(linBad) x 1 logical.

% Indexing a vector returns the orientation of the vector it indexes, not of
% the index, so every result is reshaped to a column before combining.
tf = false(numel(linBad),1);
if is_sdop
    % Content in the constant part A or in any decision row of B
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
% LIN = NONZERO_POSITIONS(PARAMS,K,IS_SDOP,NC) returns the linear positions
% within vec(C_k) that carry any content, whether through the constant part
% or through a decision variable.
% INPUTS  params, k, is_sdop, nC: as for forbidden_content.
% OUTPUT  lin: column of positions in vec(C_k), [] for no content.

if is_sdop
    % Union of the nonzeros of A and the columns of B with any nonzero
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
% PARAMS = SET_PARAM(PARAMS,K,IS_SDOP,T,NC_N,NROW_N) applies the index map T
% to parameter k. A parameter with no content is resized to the new bases,
% or left empty if that is how it was stored, since an empty parameter is
% zero.
% INPUTS  params, k, is_sdop: as for forbidden_content; T: nC_n x nC index
%         map, [] for no content; nC_n, nrow_n: numel and rows of the new
%         coefficient matrix.
% OUTPUT  params: with cell k replaced.

if is_sdop
    if isempty(T)
        % No content: a resized zero, or kept empty
        if isempty(params.A{k}),    params.A{k} = [];
        else,                       params.A{k} = sparse(nC_n,1);
        end
        if isempty(params.B{k}),    params.B{k} = [];
        else,                       params.B{k} = sparse(size(params.B{k},1),nC_n);
        end
    else
        % vec(C) = A + B'*d: A maps by T, B by T'
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
    % 'sopvar': the coefficient matrix, mapped by T on its vec
    if isempty(T)
        if isempty(params{k}),      params{k} = [];
        else,                       params{k} = sparse(nrow_n,nC_n/nrow_n);
        end
    else
        params{k} = reshape(T*sparse(params{k}(:)),nrow_n,nC_n/nrow_n);
    end
end

end
