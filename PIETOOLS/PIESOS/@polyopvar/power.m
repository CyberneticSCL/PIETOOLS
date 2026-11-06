function C = power(A,d)
% C = POWER(A,D) returns the 'polyopvar' object C representing the
% elementwise power A.^D of the distributed polynomial A.
%
% INPUTS
% - A:     m x n 'polyopvar' object;
% - d:     nonnegative integer, the order of the power;
%
% OUTPUTS
% - C:     'polyopvar' object representing A.^D, a function of a single
%          copy of the variables;
%
% NOTES
% - A.^D is the elementwise power opower(A,D,1), see
%   @polyopvar/opower.m;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - power
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
% DJ, 11/06/2026: Initial coding;

if ~isa(A,'polyopvar')
    error("The base of a power of a distributed polynomial must be a 'polyopvar' object.")
end
C = opower(A,d,1);

end
