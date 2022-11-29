%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert_PIETOOLS_PDE_terms.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PIE_out = convert_PIETOOLS_PDE(PDE,out_type)
% convert_PIETOOLS_PDE takes a PDE structure and calls the appropriate PDE
% to PIE converter function, returning the associated PIE structure.
% Input argument "type" is optional, and should be equal to ''pie'' when
% specified, as only conversion from PDE to PIE is supported.
%
% A Partial Integral Equation is defined by 12 PI operators as
%
% Tw \dot{w}(t) + Tu \dot{u}(t) + T \dot{x}(t) = A  x(t) + B1 w(t)  + B2 u(t);
%                                         z(t) = C1 x(t) + D11 w(t) + D12op u(t);
%                                         y(t) = C2 x(t) + D21 w(t) + D22op u(t);
%
% This script takes a user-defined PDE system in the format outlined in the
% header of "pde_struct" and converts it to a PIE by defining the 11
% PI operators {T,  Tw,  Tu;
%               A,  B1,  B2;
%               C1, D11, D12;
%               C2, D21, D22}; 
%
% Depending on the structure of the input "PDE", the "terms", "batch", or
% "terms_old" conversion function will be called. See the separate
% functions and the manual for more information on the input format.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, SS, DJ, 08/18/2022.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only conversion to PIE is allowed
if nargin==2 && ~strcmpi(out_type,'pie')
    error('Second argument must be ''pie'', or must be omitted. Only conversion from PDE to PIE is supported.')
end

% Call the appropriate converter function.
if isa(PDE,'pde_struct') || (isa(PDE,'struct') && isfield(PDE,'x'))
    % Convert 'struct' to 'pde_struct'.
    if isa(PDE,'struct')
        PDE = pde_struct(PDE);
    end
    % Call the "@pde_struct/convert" function.
    PDE = initialize(PDE,true);
    PIE_out = convert(PDE,'pie');
elseif isfield(PDE,'n0') || isfield(PDE,'n1') || isfield(PDE,'n2')
    PIE_out = convert_PIETOOLS_PDE_batch(PDE);
elseif isfield(PDE,'n') || isfield(PDE,'ODE')  || isfield(PDE,'PDE')
    PIE_out = convert_PIETOOLS_PDE_terms_legacy(PDE);
else
    error('The input PDE is not appropriately specified. Please define your PDE as a "pde_struct" class object, and consult the manual and examples for illustration of the structure.')
end



end