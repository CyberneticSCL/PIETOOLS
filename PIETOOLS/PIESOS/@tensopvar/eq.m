function logval = eq(A,B)
% LOGVAL = EQ(A,B) tests whether two tensor-PI operators are equal.
%
% INPUTS
% - A, B:   'tensopvar' objects of which to verify equality. Can also have
%           B=0 to test whether A is a zero operator;
%
% OUTPUTS
% - logval: logical value set to true if the operators are equal, or false
%           otherwise. Equality requires the input and output function
%           spaces of A and B to be equal, and the parameters to match. For
%           B=0, it is sufficient for only one factor in each term to be
%           zero;
%
% NOTES
% - No simplification of the operators is performed. If the terms in the
%   tensor-PI operators are the same but reordered, `eq` will not catch
%   this, and will report them to be different.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eq
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

% Assume the tensor-PI operators are not equal
logval = false;

if isa(A,'tensopvar') && isa(B,'tensopvar')
    if isequal(A,B)                                                        
        logval = true;                                                     
        return                                                             
    end
    if ~isequal(A.dims,B.dims) || ~isequal(A.vars,B.vars) ||...
            ~isequal(A.dom,B.dom) || ~isequal(A.depmat,B.depmat)
        return
    end
    C = A-B;
    logval = eq(C,0);
    return
elseif ~isa(A,'tensopvar')
    logval = eq(B,A);
    return
elseif ~isa(B,'tensopvar')
    if ~isequal(B,0)
        error("Equality can only be tested with 0 or other 'tensopvar' objects.")
    end
end

% For each term in the operator, check if at least one factor is below
% tolerance
[~,is_zeroA] = clean_op(A);
logval = is_zeroA;

end