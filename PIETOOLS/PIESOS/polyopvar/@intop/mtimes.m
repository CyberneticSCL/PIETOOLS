function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'intop' object C representing the product of
% the multiplier operator defined by A with the functional defined by B.
%
% INPUTS
% - A:      m x p 'double', 'polynomial', or 'dpvar' object (or a scalar)
%           representing some multiplier operator to multiply with B;
% - B:      p x n 'intop' object;
%
% OUTPUTS
% - C:      m x n 'intop' object representing the product of the
%           multiplier A with the functional defined by B (a 'tensopmat'
%           object if A is an 'intop' and B a 'tensopmat');
%
% NOTES
% - For a constant ('double') multiplier A, the product acts directly on
%   the sparse coefficients B.params. For a 'polynomial' or 'dpvar'
%   multiplier, the kernels of B are multiplied by A and expressed in the
%   same monomials and decision variables again.
% - If A is an 'intop' and B a 'tensopvar', C is the composition of A with
%   B (see 'mtimes_tensopvar'). If B is a 'tensopmat' with one row, C is
%   the 'tensopmat' of the compositions of A with each element of B. If B
%   is a scalar, C = B*A. Products of two 'intop' objects are not
%   supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
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

% Check that multiplication is supported
if isa(A,'intop') && isa(B,'tensopvar')
    C = mtimes_tensopvar(A,B);
    return
elseif isa(A,'intop') && isa(B,'tensopmat')
    if size(B.ops,1)>1
        error("Composition of 'intop' with multi-row 'tensopmat' is currently not supported.")
    end
    % Take the composition of A with each of the elements of B
    C = B;
    for j=1:size(B.ops,2)
        C.ops{j} = mtimes_tensopvar(A,B.ops{1,j});
    end
    C.dims.out = size(A,1);
    C.vars.out = cell(1,0);
    C.dom.out = zeros(0,2);
    C.depmat.out = zeros(1,0);
    C.chkval = true;
    return
end
if isa(A,'intop') && isa(B,'intop')
    error("Multiplication of functionals is not supported.")
elseif isa(A,'intop') && all(size(B)==1)
    % Allow right-multiplication with a scalar
    C = B*A;
    return
elseif isa(A,'intop')
    error("Left-multiplication with 'intop' objects is not supported.")
elseif ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'intop' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end
if (~B.chkval && ~isvalid(B))
    error("The functional is not properly declared.")
end

% Check that the inner dimensions of the objects match
Adim = size(A);     Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Multiplication with an identity matrix leaves the functional unchanged
if Adim(1)==Adim(2)
    if isa(A,'double')
        is_id_A = isequal(A,eye(Adim(1)));
    else
        % polynomial/dpvar multiplier (possibly depending on variables, so not an 'intop'):
        % A is the identity if A-I has no coefficients above tolerance
        D_A = A - eye(Adim(1));
        Cval = nonzeros(D_A.C);
        is_id_A = isempty(Cval) || max(abs(Cval))<=1e-12;
    end
    if is_id_A
        C = B;
        C.chkval = true;
        return
    end
end

% Multiplication with the functional amounts to multiplication with the
% parameters
if isa(A,'double')
    % Constant multiplier: only the nonzeros of the sparse parameters change
    C = B;
    nd = numel(B.Zd);
    if all(Adim==1)
        C.params = cellfun(@(a) a*A,B.params,'UniformOutput',false);
    else
        C.params = cellfun(@(a) lmult_params(A,a,nd),B.params,'UniformOutput',false);
    end
else
    % Polynomial or dpvar multiplier: multiply the kernels, and express them in
    % terms of the same monomials and decision variables again
    kernels_B = kernels(B);
    C = intop(cellfun(@(a) A*a,kernels_B,'UniformOutput',false),B.intmat,B.vars,B.dom);
end
C.chkval = true;

end
