function F = constant2polyopvar(C)
% F = CONSTANT2POLYOPVAR(C) takes a constant or decision variable C and
% converts it to a degree-0 distributed monomial F.
%
% INPUTS
% - C:  'double', 'polynomial', or 'dpvar' column vector (or 'polyopvar',
%       returned unchanged) representing a constant value or variable,
%       which cannot depend on spatial variables;
%
% OUTPUTS
% - F:  'polyopvar' object representing the same term as C, expressed as
%       a degree-0 distributed polynomial;
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
% DJ, 10/15/2026: Initial coding;

% Allow the input to already be of type 'polyopvar'
if isa(C,'polyopvar')
    F = C;
    return
end

% Determine which independent variables the term depends on
if isa(C,'double')
    pvarname = {};
    C = polynomial(C);
elseif isa(C,'polynomial') || isa(C,'dpvar')
    pvarname = C.varname;
else
    error("Input must be of type 'double', 'polynomial', or 'dpvar'.")
end
if ~isempty(pvarname)
    error("Constants depending on spatial variables are not currently supported by constant2polyopvar (domain of these variables is unknown).")
elseif size(C,2)~=1
    error("constant2polyopvar currently only supports column-vector-valued constants.")
end

% Convert coefficient to a plain (no spatial variable) multiplier operator
novars = struct('in',{cell(1,0)},'out',{cell(1,0)});
nodom = struct('in',zeros(0,2),'out',zeros(0,2));
Cop = mat2tensopvar(C,novars,nodom);

% Declare a polyopvar representation of the term, as a single degree-0
% monomial (the constant "1") with coefficient Cop
F = polyopvar();
F.varname = cell(0,1);
F.Z = sparse(1,0);
F.C = tensopmat(Cop);
F = initialize(F);

end