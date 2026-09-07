function Pop_out = uminus(Pop_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pop_out = uminus(Pop_in) or -Pop_in negates an sdopvar operator.
%
% INPUT
% Pop_in:   'sdopvar' object;
%
% OUTPUT
% Pop_out:  'sdopvar' object representing -Pop_in;
%
% NOTES:
% Negating vec(C) = A + B'*d means negating both A and B, so the variables,
% domain, dimensions, monomial bases and decision variable list are all
% unchanged. The canonical multiplier form constrains which columns may hold
% content, not the values there, so negation cannot violate it and the
% coefficients are written directly rather than rebuilt through the
% constructor.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/07/2026

Pop_out = Pop_in;
for ii=1:numel(Pop_in.params.A)
    Pop_out.params.A{ii} = -Pop_in.params.A{ii};
    Pop_out.params.B{ii} = -Pop_in.params.B{ii};
end

end
