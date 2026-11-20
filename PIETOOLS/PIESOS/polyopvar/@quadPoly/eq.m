function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two quadpoly P1 and P2 with tolerance tol
% Version: 1.0
%
% INPUT
% P1, P2: quadpoly class objects; one of them may be the scalar 0
% tol: (optional) acceptable tolerance value, default 1e-14. If every
%      coefficient of P1-P2 is smaller than tol in magnitude, then P1=P2
%
% OUTPUT
% logval: returns 1 if the objects are equal, 0 if not equal
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
% AT, 01/21/2026: Initial coding

% Default tolerance
if nargin<3
    tol=1e-14;
end    

% A scalar 0 is replaced by the zero 'quadPoly' of the size of the other
if ~isa(P1,'quadPoly')&&(P1==0)
    P1 = 0*P2;
elseif  ~isa(P2,'quadPoly')&&(P2==0) 
    P2 = 0*P1;
elseif ~isa(P1,'quadPoly')|| ~isa(P2,'quadPoly')
    error('To check equality either both values must be quadPoly objects, or one of them have to be zero');
end

diff = P1 - P2; % compute the difference between quadpoly objects

% Equal if every coefficient of the difference is below tol in magnitude
if max(max(abs(diff.C))) < tol
    logval = true;
else
    logval = false;
end

end
