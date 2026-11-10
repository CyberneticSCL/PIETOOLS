function [Cop,is_zero] = clean_op(Cop,ztol)
% [COP,IS_ZERO] = CLEAN_OP(COP,ZTOL) takes a 'tensopvar' object COP and
% removes terms (rows of COP.ops) whose tensor product is zero up to the
% specified tolerance ZTOL, i.e. terms in which at least one factor is
% zero.
%
% INPUTS
% - Cop:    'tensopvar' object representing a tensor-PI operator;
% - ztol:   (optional) scalar value specifying the tolerance below which a
%           factor's coefficients are assumed to be (collectively) zero.
%           Defaults to 1e-12;
%
% OUTPUTS
% - Cop:    'tensopvar' object representing the same tensor-PI operator as
%           the input up to the specified tolerance, with any term (row of
%           Cop.ops) in which some factor is zero up to ztol removed. At
%           least one term is always retained, even if every term is zero,
%           so that Cop.matdim/vars/dom/depmat (derived from Cop.ops)
%           remain well-defined;
% - is_zero: true if EVERY term of the input was zero up to ztol (so Cop
%           represents the zero operator and the single term retained for
%           shape is itself zero); false otherwise;
%
% NOTES
% - Each factor is cleaned via a polymorphic 'clean_op(op_ij,ztol)' call,
%   dispatched by MATLAB to @sopvar/clean_op.m, @intop/clean_op.m, or back
%   into this function, by op_ij's runtime class.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - clean_op
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

% Check that the operator is properly specified
if ~Cop.chkval
    Cop = initialize(Cop);
end

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Clean every factor (recursing into nested tensopvar/intop/sopvar),
% and identify terms (rows of Cop.ops) in which some factor is zero
ops = Cop.ops;
[ntrms,nfctrs] = size(ops);
is_zero_trm = false(ntrms,1);
for i=1:ntrms
    for j=1:nfctrs
        [ops{i,j},is_zero_ij] = clean_op(ops{i,j},ztol);
        if is_zero_ij   
            % Term's tensor product is zero regardless of the remaining
            % factors
            is_zero_trm(i) = true;
            break
        end
    end
end

% Always keep at least one term, so matdim/vars/dom/depmat stay well-defined
is_zero = ntrms>0 && all(is_zero_trm);
if is_zero
    Cop.ops = ops(1,:);
    for j=1:nfctrs
        Cop.ops{1,j} = 0*Cop.ops{1,j};
    end
else
    Cop.ops = ops(~is_zero_trm,:);
end
Cop.chkval = true;

end
