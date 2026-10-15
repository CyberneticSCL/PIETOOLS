function [f,is_zero] = clean_poly(f,ztol)
% [F,IS_ZERO] = CLEAN_POLY(F,ZTOL) takes a 'polyopvar' object F and removes
% coefficients below the specified tolerance ZTOL, monomials that do not
% contribute, and state variables that do not appear in any monomial.
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial
%           P(x) = C*Z(x);
% - ztol:   (optional) scalar value specifying the tolerance below which
%           coefficients are assumed to be zero. Defaults to 1e-12;
%
% OUTPUTS
% - f:      'polyopvar' object representing the same polynomial as the
%           input up to the specified tolerance:
%             * consecutive factors that are multipliers from R are
%               combined via combine_factors.m;
%             * duplicate monomials (rows of f.Z) are merged via
%               combine_terms.m, summing the operators acting on them --
%               this can reveal a monomial that is entirely zero
%               (duplicates that cancel), which the steps below remove;
%             * every block of f.C is cleaned via @tensopmat/clean_op.m
%               (coefficients below ztol zeroed, dead spatial variables of
%               C dropped);
%             * a monomial (row of f.Z, column of f.C.ops) that does not
%               contribute -- every block in that column is zero -- is
%               removed, unless it is the only surviving block for some
%               output row, in which case exactly one such monomial per
%               orphaned row is kept (zeroed) so that row's own shape
%               stays well-defined;
%             * trailing zero columns of f.Z are removed;
%             * a state variable (entry of f.varname) that does not
%               appear in any surviving monomial is removed from
%               varname/dims.in/depmat.in, and f.Z's indices are remapped
%               to stay contiguous;
% - is_zero: true if every coefficient of f is zero up to ztol, i.e. f
%           represents the identically-zero polynomial; false otherwise;
%
% NOTES
% - Dropping a dead monomial or state does not change the polynomial P
%   represents: apply_poly.m only reads the state values through the
%   nonzero entries of the surviving rows of Z, so a removed state is
%   never read by a surviving monomial, and a removed monomial is zero.
% - Monomial removal precedes state removal: a state may appear to be
%   used only because it appears in a monomial that is itself removed.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - clean_poly
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
% DJ, 10/15/2026: Initial coding;

if nargin==1
    ztol = 1e-12;
end

% Combine the consecutive factors that are multipliers from R, which may make monomials coincide
f = combine_factors(f);

% Merge duplicate monomials first
f = combine_terms(f);

% Clean every nonempty block of the coefficient operator, tracking which
% monomials (columns of C.ops) are now entirely zero
Cop = f.C;
ops = Cop.ops;
[nr,nZ] = size(ops);
mono_is_zero = true(1,nZ);
for j=1:nZ
    for i=1:nr
        op_ij = ops{i,j};
        if isempty(op_ij)
            continue
        end
        if isa(op_ij,'polynomial') || isa(op_ij,'dpvar') || isa(op_ij,'double')
            [ops{i,j},is_zero_ij] = clean_dpvar(op_ij,ztol);   % block acting on R, see tensopmat
        else
            [ops{i,j},is_zero_ij] = clean_op(op_ij,ztol);
        end
        mono_is_zero(j) = mono_is_zero(j) && is_zero_ij;
    end
end
Cop.ops = ops;

% Drop zero monomials, except keep exactly one per output row that would
% otherwise be left with no surviving block at all (preserves that row's
% own shape.
keep_mono = ~mono_is_zero;
for i=1:nr
    row_nonempty = ~cellfun('isempty',ops(i,:));    % isempty returns a double for a dpvar
    if any(keep_mono & row_nonempty)
        continue
    end
    j0 = find(row_nonempty,1);
    if ~isempty(j0)
        keep_mono(j0) = true;
    end
end
if ~all(keep_mono)
    Cop.ops = Cop.ops(:,keep_mono);
    Cop.dims.in = Cop.dims.in(keep_mono);
    Cop.depmat.in = Cop.depmat.in(keep_mono,:);
    f.Z = f.Z(keep_mono,:);
end
% Get rid of trailing zero columns in Z. Zeros before the last nonzero column may be the
% positions of factors that are multipliers from R, see `polyopvar`, and removing them
% would shift the factors that follow.
f.Z = f.Z(:,1:max([0,find(any(f.Z,1),1,'last')]));

% Re-clean: now that some monomials are gone, a spatial variable of C may
% have become dead too (no surviving monomial references it).
[f.C,is_zero] = clean_op(Cop,ztol);

% Drop state variables that no longer appear in any surviving monomial,
% remapping f.Z's indices to stay contiguous
Zmat = f.Z;
p = numel(f.varname);
state_used = false(p,1);
for k=1:p
    state_used(k) = nnz(Zmat==k)>0;
end
if ~all(state_used)
    remap = zeros(p,1);
    remap(state_used) = 1:nnz(state_used);
    [ridcs,cidcs,vals] = find(Zmat);
    Zmat = sparse(ridcs,cidcs,remap(vals),size(Zmat,1),size(Zmat,2));
    f.Z = Zmat;
    f.varname = reshape(f.varname(state_used),[],1);
    f.dims.in = f.dims.in(state_used);
    f.depmat.in = f.depmat.in(state_used,:);
end

f.chkval = true;

end
