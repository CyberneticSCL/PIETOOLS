function logval = eq(P1,P2,tol)
% LOGVAL = EQ(P1,P2,TOL) tests whether two 'sdopvar' objects represent the
% same decision operator, after aligning their decision-variable and
% monomial bases.
%
% INPUTS
% - P1, P2: 'sdopvar' objects of the same size; one of them may instead be
%           the scalar 0, in which case the other is tested for being zero;
% - tol:    (optional) scalar tolerance, default 1e-14. Entries differ if
%           their absolute difference is >= tol;
%
% OUTPUTS
% - logval: logical, true if P1 and P2 have the same dimensions, variables,
%           domains and number of parameters, and, once both are put on
%           common decision-variable and monomial bases, their
%           coefficients A and B agree entrywise within tol. This tests the
%           stored coefficients, which determine the operator since the
%           canonical multiplier form is a class invariant;
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
%
% Initial coding DB - 08/17/2026

% Default tolerance
if nargin<3
    tol = 1e-14;
end
% A scalar 0 operand: test the other for being the zero operator
if ~isa(P1,'sdopvar') && isequal(P1,0)
    logval = eq_zero(P2,tol);
    return
end
if ~isa(P2,'sdopvar') && isequal(P2,0)
    logval = eq_zero(P1,tol);
    return
end
if ~isa(P1,'sdopvar') || ~isa(P2,'sdopvar')
    error('Both inputs must be sdopvar objects, or one input must be zero.');
end

% Different dimensions, spaces, domains or cell counts: not equal
if any(P1.dims~=P2.dims)
    logval = false;
    return
end
if ~isequal(cellstr(string(P1.vars.in(:))),cellstr(string(P2.vars.in(:)))) || ~isequal(cellstr(string(P1.vars.out(:))),cellstr(string(P2.vars.out(:))))
    logval = false;
    return
end
if any(any(P1.dom.in~=P2.dom.in)) || any(any(P1.dom.out~=P2.dom.out))
    logval = false;
    return
end
if numel(P1.params.A)~=numel(P2.params.A)
    logval = false;
    return
end

% Put both objects in the same decision-variable and monomial coordinates.
[ops,T] = sync_basis({P1,P2});
P1 = ops{1};    P2 = ops{2};

% Compare the coefficients cell by cell on the common bases
logval = true;
for ii=1:numel(P1.params.A)
    [A1,B1] = apply_basis_map(T{1},P1.params.A{ii},P1.params.B{ii});
    [A2,B2] = apply_basis_map(T{2},P2.params.A{ii},P2.params.B{ii});
    % Test the stored nonzeros: 'all(abs(D(:))<tol)' on a sparse D would
    % build a full q*nC logical (the rows of B are the decision variables).
    logval = logval && ~any(abs(nonzeros(A1-A2))>=tol) ...
                    && ~any(abs(nonzeros(B1-B2))>=tol);
    if ~logval
        return
    end
end
end

function logval = eq_zero(P,tol)
% LOGVAL = EQ_ZERO(P,TOL) is true when every stored coefficient of the
% 'sdopvar' P, constant and decision-variable part, is below TOL in modulus.
% INPUTS  P: 'sdopvar' object (or a scalar, tested against 0); tol: scalar.
% OUTPUT  logval: logical.
if ~isa(P,'sdopvar')
    logval = isequal(P,0);
    return
end
logval = true;
for ii=1:numel(P.params.A)
    % Test the stored nonzeros only: no full q*nC logical is built
    logval = logval && ~any(abs(nonzeros(P.params.A{ii}))>=tol) ...
                    && ~any(abs(nonzeros(P.params.B{ii}))>=tol);
    if ~logval
        return
    end
end
end
