function logval = eq(P1,P2,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = eq(P1,P2,tol) tests equality of two sopvars P1 and P2 with tolerance tol
% Date: 6/13/19
% Version: 1.0
% 
% INPUT
% P1, P2: sopvar class objects
% tol: acceptable tolerance value. If max(P1-P2)<tol, then P1=P2
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
% AT, 01/21/2026: Initial coding

if nargin<3
    tol=1e-14;
end    

if ~isa(P1,'sopvar')&&(P1==0)
    P1 = 0*P2;
    % ndopvar P1
    % P1.dom = P2.dom;
    % P1.dim = P2.dim;
    % P1.vars = P2.vars;
elseif  ~isa(P2,'sopvar')&&(P2==0) 
    P2 = 0*P1;
    % ndopvar P2
    % P2.dom = P1.dom;
    % P2.dim = P1.dim;
    % P2.vars = P1.vars;
elseif ~isa(P1,'sopvar')|| ~isa(P2,'sopvar')
    error('To check equality either both values must be sopvar objects, or one of them have to be zero');
end

logval = true;
if any(P1.dims~=P2.dims)
    disp('sopvars have different dimensions and hence cannot be equal');
    logval = false;
    return
end

% Check that domain and variables match
if any(any(P1.dom_1~=P2.dom_1))|| any(any(P1.dom_2~=P2.dom_2)) || any(any(P1.dom_3~=P2.dom_3))
    disp('Operators have different intervals and hence cannot be equal');
end


 
logval = true;
for ii = 1:numel(P1.params)
    temp_log_val = eq(P1.params{ii}, P2.params{ii}, tol);
    logval = logval*temp_log_val;
end

end