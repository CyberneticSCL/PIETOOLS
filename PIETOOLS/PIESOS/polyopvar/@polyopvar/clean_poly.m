function [f,is_zero] = clean_poly(f,ztol)
% [F,IS_ZERO] = CLEAN_POLY(F,ZTOL) removes from the 'polyopvar' object F
% the coefficients below the tolerance ZTOL, the monomials that do not
% contribute, and the state variables that appear in no monomial.
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial
%           P(x) = C*Z(x);
% - ztol:   (optional) scalar, the tolerance below which coefficients are
%           set to zero. Defaults to 1e-12;
%
% OUTPUTS
% - f:      'polyopvar' object representing the same polynomial up to the
%           tolerance ztol, with consecutive multipliers from R combined,
%           duplicate monomials (rows of f.Z) merged, coefficients below
%           ztol set to zero, unused spatial variables of f.C dropped,
%           zero monomials and trailing zero columns of f.Z removed, and
%           states that appear in no monomial removed;
% - is_zero: logical, true if all coefficients of f are zero up to ztol;
%
% NOTES
% - A row-block of f.C of which all blocks are zero keeps one (zero)
%   monomial, so that its dimensions remain defined.
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
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
            [ops{i,j},is_zero_ij] = clean_dpvar(op_ij,ztol);   % block acting on R
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
% The edits above unset chkval, but C is still consistent: f was checked
Cop.chkval = true;
% Get rid of trailing zero columns in Z. Zeros before the last nonzero column may be the
% positions of factors that are multipliers from R, and removing them
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
