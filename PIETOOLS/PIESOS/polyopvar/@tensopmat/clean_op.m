function [Cop,is_zero] = clean_op(Cop,ztol)
% [COP,IS_ZERO] = CLEAN_OP(COP,ZTOL) takes a 'tensopmat' object COP and
% cleans each of its (nonempty) blocks up to the specified tolerance ZTOL.
%
% INPUTS
% - Cop:    'tensopmat' object representing a block matrix of tensor-PI
%           operators;
% - ztol:   (optional) scalar value specifying the tolerance below which a
%           block's coefficients are assumed to be (collectively) zero.
%           Defaults to 1e-12;
%
% OUTPUTS
% - Cop:    'tensopmat' object representing the same block matrix as the
%           input up to the specified tolerance, with every nonempty block
%           Cop.ops{i,j} cleaned up to ztol ('clean_dpvar' for a block of
%           type 'double', 'polynomial' or 'dpvar', 'clean_op' otherwise).
%           Any input/output variable that no block structurally depends
%           on (a column of Cop.depmat.in/.out that is entirely zero) is
%           removed from vars.in/.out, dom.in/.out, and depmat.in/.out;
% - is_zero: true if every block (already empty, or cleaned to zero) is
%           zero, i.e. Cop represents the zero operator; false otherwise;
%

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
if ~Cop.chkval && ~isvalid(Cop)
    error("Tensopmat object is not properly declared.")
end

% Set default tolerance
if nargin==1
    ztol = 1e-12;
end

% Clean every nonempty block, tracking whether every block (empty, or
% cleaned to zero) is zero
ops = Cop.ops;
[p,d] = size(ops);
is_zero = true;
for i=1:p
    for j=1:d
        op_ij = ops{i,j};
        if isempty(op_ij)
            continue   % already the zero operator for this block
        end
        if is_native_op(op_ij)
            [ops{i,j},is_zero_ij] = clean_dpvar(op_ij,ztol);
        else
            [ops{i,j},is_zero_ij] = clean_op(op_ij,ztol);
        end  
        is_zero = is_zero && is_zero_ij;
    end
end

Cop.ops = ops;

% Drop any input/output variable no block structurally depends on
rtn_vars_out = any(Cop.depmat.out,1);
if ~all(rtn_vars_out)
    Cop.vars.out = Cop.vars.out(rtn_vars_out);
    Cop.dom.out = Cop.dom.out(rtn_vars_out,:);
    Cop.depmat.out = Cop.depmat.out(:,rtn_vars_out);
end
rtn_vars_in = any(Cop.depmat.in,1);
if ~all(rtn_vars_in)
    Cop.vars.in = Cop.vars.in(rtn_vars_in);
    Cop.dom.in = Cop.dom.in(rtn_vars_in,:);
    Cop.depmat.in = Cop.depmat.in(:,rtn_vars_in);
end

Cop.chkval = true;

end
