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
%               Alternatively, 'char' or cell of 'char' objects with the
%               flags 'silent' (suppress the display), 'nonlinear' (return
%               an 'npie_struct' also for a linear PDE) and/or 'Top' of
%               "convert_PIETOOLS_PDE", which are passed on to it;
%
% OUTPUT:
% - PIE:        A pie_struct object (or an npie_struct object, if the PDE
%               is nonlinear or the flag 'nonlinear' is specified) defining
%               the equivalent PIE representation of the input PDE, if this
%               representation exists.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pde_struct/convert
%
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
% Initial coding DJ - 10/11/2022
% Call combines 1D-2D converter, DJ - 10/17/2024
% Display successful connversion message DB -12/22/24
% DJ, 06/01/2025: Move pre-processing of PDE to "convert_PIETOOLS_PDE";
% DJ, 11/02/2026: Accept the flags 'silent', 'nonlinear', 'Top' of
%                   "convert_PIETOOLS_PDE" as second input;

% We support only PDE to PIE conversion. The second argument is 'pie', or
% flags of convert_PIETOOLS_PDE, with the 'pie' left out.
if nargin<2
    out_type = {};
end
flag = cellstr(out_type);
if ~all(ismember(flag,{'pie','silent','nonlinear','Top','return_Top'}))     % DJ, 11/02/2026
    error(['Only conversion from PDE to PIE is supported. The second argument must be ''pie'', ',...
        'or flags of convert_PIETOOLS_PDE (''silent'', ''nonlinear'', ''Top''), or be omitted.'])
end
flag = flag(~strcmp(flag,'pie'));

% Call the PDE to PIE converter.
PIE = convert_PIETOOLS_PDE(PDE,flag);

% DB 12/22/24
if ~ismember('silent',flag)                                                 % DJ, 11/02/2026
    fprintf('\n --- Conversion to PIE was successful --- \n');
end
end