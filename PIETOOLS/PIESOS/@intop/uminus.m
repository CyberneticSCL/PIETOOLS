function C_out = uminus(C_in)
% C_OUT = UMINUS(C_IN) returns the 'intop' object C_OUT representing -C,
% for C the functional represented by the 'intop' object C_IN.
%
% INPUTS
% - C_in:   'intop' object;
%
% OUTPUTS
% - C_out:  'intop' object representing the negative of the input
%           functional;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - uminus
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

% Check that the functional is properly specified
if ~C_in.chkval && ~isvalid(C_in)
    error("The functional is not properly declared.")
end

% The kernels defining the negative functional are just the negative
% kernels
C_out = C_in;
C_out.params = cellfun(@(a)-a,C_in.params,'UniformOutput',false);
C_out.chkval = true;

end