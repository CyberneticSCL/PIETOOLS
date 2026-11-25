function logval = eq(A,B)
% LOGVAL = EQ(A,B) tests whether two tensor-PI operators are equal.
%
% INPUTS
% - A, B:   'tensopvar' objects of which to verify equality. Can also have
%           B=0 (or A=0) to test whether A (B) is a zero operator;
%
% OUTPUTS
% - logval: 'logical' scalar, true if the operators are equal, and false
%           otherwise. Equality requires the input and output function
%           spaces of A and B to be equal, and A-B to be zero. For B=0, it
%           is sufficient for one factor in each term to be zero (up to
%           1e-12);
%
% NOTES
% - No simplification of the operators is performed. Single-factor
%   operators are compared by their parameters, but operators with
%   several factors are only found equal if they are identical, or if
%   both are zero; e.g. the same terms in a different order are reported
%   to be different.
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 10/15/2026: Initial coding;

% Assume the tensor-PI operators are not equal
logval = false;

% Two operators: equal if identical, or if the function spaces match
% and A-B is zero
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
% Otherwise the operator is compared with 0
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