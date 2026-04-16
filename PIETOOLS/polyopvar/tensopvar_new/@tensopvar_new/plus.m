function C = plus(A,B)
% C = PLUS(A,B) returns the 'tensopvar' object C representing the sum
% of the tensor-PI operators defined by the 'tensopvar' objects A and B
%
% INPUTS
% - A:      m x n 'tensopvar' object
% - B:      m x n 'tensopvar' object
%
% OUTPUTS
% - C:      m x n 'tensopvar' object representing the sum of the
%           tensor-PI operators defined by A and B
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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
% DJ, 04/15/2026: Initial coding

% Check that the inputs are appropriate
if ~isa(A,'tensopvar_new')
    if ~isa(A,'nopvar') && ~isa(A,'ndopvar')
        error("Addition of 'tensopvar' objects with non-'tensopvar' objects is not supported.")
    else
        A = ndopvar2tensopvar_new(A);
    end
elseif ~isa(B,'tensopvar_new')
    if ~isa(B,'nopvar') && ~isa(B,'ndopvar')
        error("Addition of 'tensopvar' objects with non-'tensopvar' objects is not supported.")
    else
        B = ndopvar2tensopvar_new(B);
    end
end

% Make sure the matrix dimensions of the operators match
if ~isequal(A.dims,B.dims) || any(A.type~=B.type)
    error("Matrix dimensions of each factor in the tensopvar objects must match.")
end
% Make sure the input and output variables match
if ~isequal(pvar2varname(A.vars),pvar2varname(B.vars)) || ~isequal(A.dom,B.dom)
    error("Spatial variables and domains of operators must match.")
end
% Make sure the individual operators depend on the same variables
if ~isequal(A.depmat1,B.depmat1) || ~isequal(A.depmat2,B.depmat2)
    error("Spatial variables in each factor in the tensopvar objects must match.")
end

% The operators defining the sum are given by the sum of the
% operators
C = A;
if isscalar(A.ops)
    % The tensor-PI operator is really just a PI operator
    % --> we can use the PI addition routine
    C.ops{1} = A.ops{1}+B.ops{1};
else
    % The tensor-PI operator is defined by multiple factors
    % --> we have no closed-form expression for the sum (since the sum of 
    %       the product does not equal product of sum)
    % Instead, we concatenate the operators to represent their sum
    C.ops = [A.ops; B.ops];
end

end