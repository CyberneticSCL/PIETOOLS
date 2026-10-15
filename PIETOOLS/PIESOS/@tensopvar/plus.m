function C = plus(A,B)
% C = PLUS(A,B) returns the 'tensopvar' object C representing the sum
% of the tensor-PI operators defined by the 'tensopvar' objects A and B.
%
% INPUTS
% - A:      m x n 'tensopvar' object (or 0, or empty);
% - B:      m x n 'tensopvar' object (or 0, or empty);
%
% OUTPUTS
% - C:      m x n 'tensopvar' object representing the sum of the
%           tensor-PI operators defined by A and B;
%
% NOTES
% - If A or B is 0 or empty, the other is returned.
% - If A and B both have a single factor, the factors are added.
% - If A and B have factors of the same dimensions and dependencies, their
%   terms (rows of ops) are concatenated.
% - Otherwise, A and B are kept as the two terms of a nested operator, with
%   a single factor {A;B}.

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
% DJ, 10/15/2026: Initial coding;

% Check that the inputs are appropriate
if isequal(A,0) || isempty(A)
    C = B;
    return
elseif isequal(B,0) || isempty(B)
    C = A;
    return
end
if ~isa(A,'tensopvar')
    try A = tensopvar(A);
    catch
        error("Object could not be converted to 'tensopvar'; addition is not supported.")
    end
end
if ~isa(B,'tensopvar')
    try B = tensopvar(B);
    catch
        error("Object could not be converted to 'tensopvar'; addition is not supported.")
    end
end

% Make sure the objects are properly declared
if ~A.chkval
    A = initialize(A);
end
if ~B.chkval
    B = initialize(B);
end

% Make sure the matrix dimensions of the operators match
if ~isequal(size(A),size(B))
    error("Matrix dimensions of each factor in the tensopvar objects must match.")
end
% Make sure the input and output variables match
if ~isequal(A.vars,B.vars) || ~isequal(A.dom,B.dom)
    error("Spatial variables and domains of operators must match.")
end
% Make sure the degrees in the input and output variables match
if ~isequal(A.deps,B.deps)
    error("Input and output domains of the operators must match.")
end

% Distinguish three cases
% Single-factor branch only if B too has a single factor: B.ops{i,1} is
% all that is added, so any further factors of B would be dropped
if size(A.ops,2)==1 && size(B.ops,2)==1
    % The tensor-PI operators consist of only a single factor
    % --> we can just add the two factors
    C = A;
    C.ops = A.ops(1);
    for i=2:size(A.ops,1)
        C.ops{1} = C.ops{1} + A.ops{i,1};
    end
    for i=1:size(B.ops,1)
        C.ops{1} = C.ops{1} + B.ops{i,1};
    end
elseif isequal(A.dims,B.dims) && isequal(A.depmat,B.depmat) && isequal(A.otype,B.otype)
    % The dimensions and domains of the factors in A match those in B
    % --> we can concatenate the factors to represent their sum
    C = A;
    C.ops = [A.ops; B.ops];
else
    % The different operators have different structures
    % --> retain their structure within the new operator
    C = tensopvar();
    C.ops = {A;B};
    C.otype = [2,2];
    C = initialize(C);
end

end
