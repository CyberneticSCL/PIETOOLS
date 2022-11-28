function PIE = convert(PDE,out_type)
% PIE = convert(PDE,out_type)
% Convert a 'pde_struct' object PDE to an associated 'pie_struct' object
% PIE.
%
% INPUT:
% - PDE:        A pde_struct object defining a PDE in the terms format
%               (see also "@pde_struct/initialize").
% - out_type:   char object specifying to what type of system the PDE is to
%               be converted. Only out_type='pie' is supported, as we only
%               convert PDEs to PIEs. This second argument is only included
%               to match "convert" function for 'sys' type objects.
%
% OUTPUT:
% - PIE:        A pie_struct object defining the equivalent PIE 
%               representation of the input PDE, if this representation
%               exists.
%
% NOTES:
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or D. Jagt at djagt@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pde_struct/convert
%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 10/11/2022

% We support only PDE to PIE conversion.
if nargin==2 && ~strcmpi(out_type,'pie')
    error('Only conversion from PDE to PIE is supported. The second argument must be set to ''pie'', or be omitted.')
end

% % % Get rid of any higher-order temporal derivatives and delays in the
% % % system, and try to reduce the number of spatial variables.
PDE = initialize(PDE,true);
if PDE.has_hotd
    fprintf(['\n','Higher-order temporal derivatives were encounterd:\n']);
    fprintf([' --- Running "expand_tderivatives" to reduce to first-order temporal derivatives ---\n']);
    PDE = expand_tderivatives(PDE);
end
if PDE.has_delay
    fprintf(['\n','Delayed states on inputs were encounterd:\n']);
    fprintf([' --- Running "expand_delays" to remove the delays --- \n']);
    PDE = expand_delays(PDE);
end
if PDE.dim>2
    fprintf(['\n','The PDE involves more than 2 spatial variables:\n']);
    fprintf([' --- Attempting to reduce the dimensionality using "combine_vars" ---\n']);
    PDE = combine_vars(PDE);
    if PDE.dim>2
        error('The number of spatial variables could not be reduced to 2 or less; conversion to PIE is not currently supported.')
    end
end

% Reorder state components, inputs, and outputs, so that these may be
% represented using PI operators.
fprintf(['\n',' --- Reordering the state components to allow for representation as PIE ---\n']);
[PDE] = reorder_comps(PDE,{'x','u','w','y','z'});

% Check the spatial dimension of the system, and call the corresponding 
% converter.
if PDE.dim<=1
    PIE = convert_PIETOOLS_PDE_terms(PDE);    
elseif PDE.dim==2
    PIE = convert_PIETOOLS_PDE_2D(PDE);
else
    error('PIETOOLS does not currently support conversion of PDEs in more than 2 spatial variables'); 
end
if isa(PIE,'struct')
    PIE = pie_struct(PIE);
end

end