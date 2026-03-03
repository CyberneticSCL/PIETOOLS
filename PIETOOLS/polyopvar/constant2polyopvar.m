function F = constant2polyopvar(C)
% F = CONSTANT2POLYOPVAR(C) takes a constant or decision variable C and
% converts it to a degree-0 distributed monomial F
%
% INPUTS
% - C:  'double', 'polynomial', or 'dpvar' object representing a constant
%       value or non-distributed polynomial function or variable;
%
% OUTPUS
% - F:  'polyopvar' object representing the same term as C, but now
%       expressed as a degree-0 distributed polynomial;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - constant2polyopvar
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
% DJ, 03/01/2026: Initial coding

% Allow the input to already be of type 'polyopvar'
if isa(C,'polyopvar')
    F = C;
    return
end

% Determine which independent variables the term depends on
if isa(C,'double')
    pvarname = {};
elseif isa(C,'polynomial') || isa(C,'dpvar')
    pvarname = C.varname;
else
    error("Input must be of type 'double', 'polynomial', or 'dpvar'.")
end

% Declare a polyopvar representation of the term
F = polyopvar();
F.pvarname = pvarname;
F.varmat = false(0,numel(pvarname));
F.degmat = zeros(1,numel(pvarname));
F.C.ops{1} = C;




end