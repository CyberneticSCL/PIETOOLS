function [vars,dom,ZL,ZR,rowIdx,colIdx,changed] = canonical_var_order(vars,dom,ZL,ZR,dims)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [vars,dom,ZL,ZR,rowIdx,colIdx,changed] = ...
%                        canonical_var_order(vars,dom,ZL,ZR,dims)
% puts the spatial variables of an sopvar/sdopvar into the canonical order
%
%       vars.out = [S2, S3]        vars.in = [S3, S1]
%
% where S1 = vars.in \ vars.out are the variables integrated out, S2 =
% vars.out \ vars.in those introduced, and S3 the pass-through variables
% common to both. S1, S2 and S3 each keep the sorted order that 'setdiff'
% and 'intersect' produce.
%
% INPUT
% vars: struct with fields 'in' and 'out', cellstr of variable names;
% dom:  struct with fields 'in' and 'out', one row per corresponding variable;
% ZL:   1 x numel(vars.out) cell of exponent columns;
% ZR:   1 x numel(vars.in)  cell of exponent columns;
% dims: 1 x 2 [rows, columns], the numbers of output and input components;
%
% OUTPUT
% vars,dom,ZL,ZR: the same content, reordered into the canonical order;
% rowIdx,colIdx:  gather indices for the ROW and COLUMN index of the
%                 coefficient matrix, so that
%                     C_new = C_old(rowIdx,colIdx),
%                 or in vectorized form
%                     vec_new(p) = vec_old(rowIdx(r) + (colIdx(c)-1)*p_rows).
%                 Both are EMPTY when the input was already canonical, so a
%                 caller can skip the gather rather than perform an identity;
% changed:        false when the input was already canonical;
%
% NOTES
% ZL{i} is bound positionally to vars.out(i) and ZR{i} to vars.in(i), so
% reordering the variables must reorder those cells and permute the monomial
% index of every stored coefficient to match. The monomial vector is
% kron(Z{1},...,Z{N}) with the FIRST variable outermost, so a permutation of
% the cells is a tensor transposition of that index.
%
% The permutation is obtained by building the full degree table both ways and
% matching rows, rather than by stride arithmetic. The table has one row per
% monomial and one column per direction, i.e. it lives on the small axis,
% and getting a tensor stride wrong here would be silent.
%
% The parameter CELL array is NOT permuted and must not be. Its multi-index
% runs over S3 in the order 'intersect(vars.in,vars.out)' returns, which this
% reordering leaves alone; a direction is recovered from a multi-index
% component with ismember against the stored order at the point of use, so
% it follows the new order.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
% Initial coding MMP, 09/09/2026

% Variable lists as rows
vin  = vars.in(:).';        vout = vars.out(:).';

% Normalize the two input shapes the class tolerates elsewhere (row
% exponent vectors, [] for an empty domain side): the code below reads the
% column form.
ZL = reshape(ZL,1,[]);      ZR = reshape(ZR,1,[]);
for i = 1:numel(ZL),    ZL{i} = ZL{i}(:);    end
for i = 1:numel(ZR),    ZR{i} = ZR{i}(:);    end
if isempty(dom.in),     dom.in  = zeros(0,2);   end
if isempty(dom.out),    dom.out = zeros(0,2);   end

% A repeated name would be silently swallowed: setdiff and intersect unique()
% their arguments, so the target lists come back shorter than the input and
% the gather below truncates the operator, discarding the coefficients of
% every direction past the duplicate.
if numel(unique(vin))~=numel(vin) || numel(unique(vout))~=numel(vout)
    error("canonical_var_order: a spatial variable name is repeated within "...
         +"vars.in or vars.out; each names one variable.")
end
% Likewise a basis cell of the wrong length would truncate rather than error,
% since ZL{i} is bound positionally to vars.out(i).
if numel(ZL)~=numel(vout) || numel(ZR)~=numel(vin)
    error("canonical_var_order: ZL and ZR must have one entry per variable "...
         +"in vars.out and vars.in (got %d/%d and %d/%d).", ...
          numel(ZL),numel(vout),numel(ZR),numel(vin))
end

% The partition: S1 input only, S2 output only, S3 shared, each sorted
S1 = reshape(setdiff(vin,vout),1,[]);
S2 = reshape(setdiff(vout,vin),1,[]);
S3 = reshape(intersect(vin,vout),1,[]);

% A pass-through variable is ONE variable, so its interval must agree on the
% two sides; operations on operators compare domains BETWEEN operands, never
% within one object, so the check is made here. Checked before the
% fast-path return below, so an already-canonical object is validated too.
[~,pl] = ismember(S3,vout);
[~,pr] = ismember(S3,vin);
% NaN~=NaN: two NaNs in the same position agree, as for isequaln, so an
% unbounded or placeholder domain passed to both sides is accepted.
dl = dom.out(pl,:);     dr = dom.in(pr,:);
bad = find(any(dl~=dr & ~(isnan(dl)&isnan(dr)),2),1);
if ~isempty(bad)
    error("canonical_var_order: variable '"+string(S3{bad})+"' appears in "...
         +"both the input and output space but is given the domain "...
         +mat2str(dom.in(pr(bad),:))+" on the input side and "...
         +mat2str(dom.out(pl(bad),:))+" on the output side.")
end

% The canonical orders
outTarget = [S2,S3];
inTarget  = [S3,S1];

changed = ~isequal(vout,outTarget) || ~isequal(vin,inTarget);
if ~changed
    % Already canonical, the common case (S1 and S2 empty, the variables in
    % sorted order): no gather, so this path stays cheap.
    rowIdx = [];    colIdx = [];
    return
end

% New position k holds the variable that was at old position ord(k).
[okL,ordL] = ismember(outTarget,vout);
[okR,ordR] = ismember(inTarget,vin);
if ~all(okL) || ~all(okR)
    error("canonical_var_order: the variable lists are inconsistent.")
end

% Monomial permutations of the reordered bases, then the reordered data
pL = monomial_gather(ZL,ordL);
pR = monomial_gather(ZR,ordR);

vars.out = outTarget;               vars.in  = inTarget;
dom.out  = dom.out(ordL,:);         dom.in   = dom.in(ordR,:);
ZL = ZL(ordL);                      ZR = ZR(ordR);

% Lift each monomial permutation over the matrix dimension. The row index of
% the coefficient is (matrix row outer, ZL monomial inner) and the column
% index likewise with ZR, so row (i,z) takes old row (i,pL(z)).
NL = prod([cellfun(@numel,ZL),1]);
NR = prod([cellfun(@numel,ZR),1]);
rowIdx = reshape(pL(:) + (0:dims(1)-1)*NL,[],1);
colIdx = reshape(pR(:) + (0:dims(2)-1)*NR,[],1);

end


