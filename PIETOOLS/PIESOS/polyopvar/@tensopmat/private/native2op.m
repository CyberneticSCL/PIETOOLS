function Top = native2op(M,X,i)
% TOP = NATIVE2OP(M,X,I) converts the matrix M of type 'double' or
% 'polynomial', which is a block of row I of the 'tensopmat' X, to the
% multiplier operator from R to the output space of row I of X, a
% 'tensopvar' object, so that the product with another 'tensopvar' (otimes,
% mtimes) is taken by the regular product of 'tensopvar' objects. The
% output space is that of the row, depmat.out(i,:), also if M depends on
% fewer of its variables (see `tensopmat`).
%
% INPUTS
% - M:      m x n matrix of type 'double' or 'polynomial';
% - X:      'tensopmat' object, with M a block in its row I;
% - i:      scalar integer, index of the row-block of X containing M;
%
% OUTPUTS
% - Top:    m x n 'tensopvar' object representing the multiplier operator
%           from R^n to the output space of row I of X;
%
% NOTES
% - A matrix of type 'dpvar' cannot be converted, as tensor-PI operators
%   cannot depend on decision variables.
% - The row must not map to a tensor product of copies of a variable
%   (X.depmat.out(i,:) <= 1).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - native2op
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

if isa(M,'dpvar')
    error("A block of type 'dpvar' cannot be converted to an operator, as tensor-PI operators cannot depend on decision variables.")
end
dpr = X.depmat.out(i,:);
if any(dpr>1)
    error("A block of type 'double' or 'polynomial' cannot be converted to an operator if its row maps to a tensor product of copies of a variable.")
end
is_out = logical(dpr);
vars = struct('in',{cell(1,0)},'out',{X.vars.out(is_out)});
dom = struct('in',zeros(0,2),'out',X.dom.out(is_out,:));
Top = mat2tensopvar(M,vars,dom);

end
