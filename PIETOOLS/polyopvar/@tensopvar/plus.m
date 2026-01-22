function C = plus(A,B)
% C = PLUS(A,B) returns the 'tensopvar' object C representing the sum
% of the coefficient operators defined by the 'tensopvar' objects A and B
%
% INPUTS
% - A:      m x n 'tensopvar' object
% - B:      m x n 'tensopvar' object
%
% OUTPUS
% - C:      m x n 'tensopvar' object representing the sum of the
%           coefficient operators defined by A and B
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
% DJ, 01/16/2026: Initial coding

% Make sure the size of the coefficients match
if any(size(A.ops)~=size(B.ops))
    error("Coefficient operators must be of same dimension for addition.")
end

% The operators defining the sum are given by the sum of the
% operators
C = A;
for ii=1:numel(A.ops)
    if isempty(A.ops{ii})
        C.ops{ii} = B.ops{ii};
    elseif isempty(B.ops{ii})
        C.ops{ii} = A.ops{ii};
    elseif isa(A.ops{ii},'nopvar') || isa(A.ops{ii},'ndopvar')
        % If the monomial is just linear, Z_{i}(x) = xk for some k,
        % then we can use the 'nopvar' plus routine to compute the sum of
        % the operators
        %   Ai*xk + Bi*xk = (Ai+Bi)*xk;
        C.ops{ii} = A.ops{ii} + B.ops{ii};
    else
        % If the monomial is more complicated, e.g. Z_{i}(x) = xk*yl,
        % we cannot simply add the coefficient
        %   (Ai{1}*xk)(Ai{2}*yk) + (Bi{1}*xk)(Bi{2}*yk)
        C.ops{ii} = [A.ops{ii}; B.ops{ii}];
    end
end

end