function logval = is_eye(A,n,ztol)
% LOGVAL = IS_EYE(A,N,ZTOL) checks whether A is an identity operator: an
% identity matrix of type 'double', 'polynomial' or 'dpvar', the
% 'polyopvar' object that maps the states to themselves, or an identity
% operator of a class with a method 'is_id_op', e.g. 'tensopvar',
% 'tensopmat' or 'intop'.
%
% INPUTS
% - A:      m x m array of type 'double' (or 'logical'), 'polynomial' or
%           'dpvar', or a 'polyopvar' object, or an operator, e.g. a
%           'tensopvar', 'tensopmat' or 'intop' object;
% - n:      (optional) nonnegative integer, the expected size of the
%           identity operator: size(A,1) = size(A,2) = n, or size(A,1) = n
%           for a 'polyopvar' object. Defaults to [], in which case any
%           size is accepted;
% - ztol:   (optional) nonnegative scalar tolerance below which
%           coefficients are taken to be zero. Defaults to 1e-14;
%
% OUTPUTS
% - logval: boolean specifying whether (true) or not (false) A is an
%           identity operator, of size n if n is specified. A scalar 1 is
%           the identity operator of any size;
%
% NOTES
% - An array is the identity if A-eye(m) is zero up to ztol.
% - A 'polyopvar' object is the identity if it is the vector of all its
%   states, in the order of A.varname, so that it represents the map x -> x
%   (see 'ispvar'). An adjoint is never the identity.
% - For an operator, the size n is checked by this function, and whether it
%   is the identity by the method 'is_id_op' of its class.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_eye
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
% DJ, 11/07/2026: Initial coding;

% Check the expected size
if nargin<2
    n = [];
elseif ~isempty(n) && (~isnumeric(n) || ~isscalar(n) || n<0 || n~=round(n))
    error("The expected size of the identity operator must be specified as a nonnegative integer.")
end
% Check the tolerance
if nargin<3
    ztol = 1e-14;
elseif ~isempty(ztol) && (~isnumeric(ztol) || ~isscalar(ztol) || ztol<0)
    error("The tolerance must be specified as a nonnegative scalar.")
end

% Assume A is not an identity operator until we have verified
logval = false;

if isnumeric(A) || islogical(A) || isa(A,'polynomial') || isa(A,'dpvar')
    % Only a square array can be the identity, of size n unless A is scalar
    m = size(A,1);
    if size(A,2)~=m || m==0 || (m>1 && ~isempty(n) && m~=n)
        return
    end
    % A-eye(m) is zero up to ztol if its coefficients are
    if isempty(ztol)
        ztol = 1e-14;   
    end
    D = A - eye(m);
    if isa(D,'dpvar')
        coeffs = nonzeros(D.C);
    elseif isa(D,'polynomial')
        coeffs = nonzeros(D.coeff);
    else
        coeffs = nonzeros(D);
    end
    logval = isempty(coeffs) || max(abs(coeffs))<=ztol;
elseif isa(A,'polyopvar')
    % The map x -> x: the vector of all states, in order
    if A.is_adjoint || isempty(A.C.ops) || (~isempty(n) && size(A,1)~=n)
        return
    end
    [is_state,xnames] = ispvar(A);
    logval = is_state && isequal(xnames(:),A.varname(:));
elseif ismethod(A,'is_id_op')
    % An operator: check the size here, and the identity by its class method
    if ~isempty(n) && ~all(size(A)==1) && ~(size(A,1)==n && size(A,2)==n)
        return
    end
    if isempty(ztol)
        logval = is_id_op(A);
    else
        logval = is_id_op(A,ztol);
    end
else
    error("The identity operator can only be checked for an object of type 'double', 'polynomial', 'dpvar' or 'polyopvar', or of a class with a method 'is_id_op'.")
end

end
