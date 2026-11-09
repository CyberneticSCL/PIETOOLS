function PIE = initialize_PIETOOLS_PIE(PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = initialize_PIETOOLS_PIE(PIE) takes a PIE data structure and
% checks that all the necessary fields are appropriately specified, and
% assigns a default value to all the optional fields that have not been
% specified.
% 
% INPUT
%   PIE: "struct", "pie_struct" or "npie_struct" class object. A "struct"
%        should have the operators of the PIE as fields, T and A
%        (mandatory) and any of Tw, Tu, B1, B2, C1, C2, D11, D12, D21, D22
%        (or Bw, Bu, Cz, Cy, Dzw, Dzu, Dyw, Dyu, K, Kw, Ku), as 'copvar',
%        'opvar' or 'opvar2d' objects. Fields vars and dom are ignored, as
%        the variables and domain are determined by T. If the struct has a
%        field f, gz, or gy, it represents a nonlinear PIE, with
%        'polyopvar' objects f (mandatory), gz and gy instead of A, Cz and
%        Cy, the other operators as 'tensopmat' objects (or e.g. 'opvar'),
%        and optionally the field space, e.g. to name the state components.
%
% OUTPUT
%   PIE: initialized "pie_struct" class object describing the same system
%        as the input, with the operators as 'copvar' objects, or, for a
%        nonlinear PIE, an initialized "npie_struct" object.
%
% NOTES:
% The initialization is done by the "initialize" method of "pie_struct",
% which throws an error if the spaces of an operator do not match those of 
% the others.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding SS - 08/01/2022
% DJ, 12/07/2024: Use new default vars s1 and s1_dum;
% DJ, 11/02/2026: Update to new pie_struct class, and support npie_struct;

if nargin==0
    PIE = pie_struct();
    return
end
if isa(PIE,'pie_struct') || isa(PIE,'npie_struct')
    PIE = initialize(PIE);
    return
elseif ~isa(PIE,'struct')
    error('Input must be a ''struct'', ''pie_struct'' or ''npie_struct'' class object.')
end

% Collect the fields of a struct in the pie_struct or npie_struct object,
% and call the internal initialize.
ops = {'T','Tw','Tu','A','B1','B2','C1','C2','D11','D12','D21','D22','Bw','Bu','Cz','Cy','Dzw','Dzu','Dyw','Dyu','K','Kw','Ku'};
PIE_in = PIE;      PIE = pie_struct();
if any(isfield(PIE_in,{'f','gz','gy'}))
    % Nonlinear PIE: the polynomials f, gz, gy replace the operators A, Cz, Cy.
    ops = {'T','Tw','Tu','Bw','Bu','Dzw','Dzu','Dyw','Dyu','B1','B2','D11','D12','D21','D22','K','Kw','Ku','f','gz','gy','space'};
    PIE = npie_struct();
end
for fname = ops(isfield(PIE_in,ops))
    PIE.(fname{1}) = PIE_in.(fname{1});
end
PIE = initialize(PIE);

end
