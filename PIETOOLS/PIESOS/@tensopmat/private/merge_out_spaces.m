function [A,B] = merge_out_spaces(A,B)
% [A,B] = MERGE_OUT_SPACES(A,B) takes two 'tensopmat' objects A, B
% expressed in terms of the same variables (see `common_vars`), and sets
% the output dependency arrays A.depmat.out(i,:), B.depmat.out(i,:) of each
% row i that consists of blocks of type 'double', 'polynomial' or 'dpvar'
% (and empty blocks) in A or B, so that A and B can be added or
% concatenated. Such a row maps to the variables appearing in its blocks,
% see `tensopmat`, and the same function can be expressed on more variables
% than appear in it:
% - If the row is of this type in both A and B, then the output variables
%   are the union of those of A and B;
% - If it is in only one of A and B, the output variables are those of the
%   other, provided these include the variables of the first.
% Other rows are not changed, nor are rows for which these conditions fail,
% for which the caller finds that the output spaces do not match.
%
% INPUTS
% - A, B:   'tensopmat' objects expressed in terms of the same variables,
%           with the same number of row-blocks. If the number of row-blocks
%           differs, or the dependency arrays A.depmat.out and B.depmat.out
%           already agree, nothing is changed;
%
% OUTPUTS
% - A, B:   'tensopmat' objects representing the same operators as the
%           inputs, with A.depmat.out(i,:) = B.depmat.out(i,:) for each row
%           i that is merged;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - merge_out_spaces
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

if size(A.ops,1)~=size(B.ops,1) || isequal(A.depmat.out,B.depmat.out)
    return      % nothing to merge if the output spaces already agree (the usual case)
end
is_nat_A = native_mask(A.ops);        is_nat_B = native_mask(B.ops);
is_empty_A = cellfun('isempty',A.ops);           is_empty_B = cellfun('isempty',B.ops);
isrow_A = any(is_nat_A,2) & all(is_nat_A | is_empty_A,2);
isrow_B = any(is_nat_B,2) & all(is_nat_B | is_empty_B,2);
for i=find(isrow_A | isrow_B)'
    a = A.depmat.out(i,:);      b = B.depmat.out(i,:);
    if isrow_A(i) && isrow_B(i)
        m = max(a,b);
    elseif isrow_A(i) && all(a<=b)
        m = b;
    elseif isrow_B(i) && all(b<=a)
        m = a;
    else
        continue
    end
    A.depmat.out(i,:) = m;      B.depmat.out(i,:) = m;
end

end
