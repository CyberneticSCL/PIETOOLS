function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOGVAL = EQ(P1,P2,TOL) tests whether two 'sopvar' objects represent the
% same operator: P1 == P2 is true when every coefficient of P1 - P2 is
% smaller than TOL in modulus.
%
% INPUT
% - P1, P2: 'sopvar' objects; one of them may be the scalar 0, which stands
%           for the zero operator on the spaces of the other;
% - tol:    (optional) scalar, the absolute tolerance on the coefficients
%           of the difference. Defaults to 1e-14;
%
% OUTPUT
% - logval: logical, true if P1 and P2 represent the same operator to the
%           tolerance, false otherwise (also when the dimensions differ).
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
% AT, 01/21/2026: Initial coding
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): help block in the
%                  house format, step comments, stale comments removed
%                  (commented-out field assignments of another class);
%                  code unchanged.

% Default tolerance
if nargin<3
    tol=1e-14;
end

% A scalar 0 stands for the zero operator on the spaces of the other operand
if ~isa(P1,'sopvar')&&(P1==0)
    P1 = 0*P2;
elseif  ~isa(P2,'sopvar')&&(P2==0)
    P2 = 0*P1;
elseif ~isa(P1,'sopvar')|| ~isa(P2,'sopvar')
    error('To check equality either both values must be sopvar objects, or one of them have to be zero');
end

% Operators of different dimensions are not equal
logval = true;
if any(P1.dims~=P2.dims)
    disp('sopvars have different dimensions and hence cannot be equal');
    logval = false;
    return
end

% Report differing domains (the difference below then fails in 'plus')
if any(any(P1.dom_1~=P2.dom_1))|| any(any(P1.dom_2~=P2.dom_2)) || any(any(P1.dom_3~=P2.dom_3))
    disp('Operators have different intervals and hence cannot be equal');
end

% Compare through the difference, so that 'plus' aligns the monomial bases: % MMP, 08/29/2026
% two objects can represent one operator on different bases ('ctranspose'   % MMP, 08/29/2026
% enlarges them where a multiplier and an integral direction share a basis). % MMP, 08/29/2026
D = P1 + uminus(P2);                                                        % MMP, 08/29/2026

% Equal iff every coefficient of the difference is below the tolerance
logval = true;
for ii = 1:numel(D.params)                                                  % MMP, 08/29/2026
    dv = nonzeros(D.params{ii});    % 'max(...,''all'')' rejects sparse input % MMP, 08/29/2026
    if ~isempty(dv) && max(abs(dv))>=tol                                    % MMP, 08/29/2026
        logval = false;                                                     % MMP, 08/29/2026
        return                                                              % MMP, 08/29/2026
    end                                                                     % MMP, 08/29/2026
end

end