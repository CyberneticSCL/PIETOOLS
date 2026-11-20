%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize_PIETOOLS_PDE.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PDE_out = initialize_PIETOOLS_PDE(PDE,silent_initialize)
% Initializes and checks the PDE formulation.
% Depending on the structure of the input "PDE", the "terms", "batch", or
% "terms_old" initialization function will be called.
%
% INPUT
% - PDE:    'pde_struct' object, or 'struct' in the terms format (field
%           x), the batch format (fields n0, n1, n2) or the legacy terms
%           format (fields n, ODE, PDE); other inputs throw an error;
% - silent_initialize:  (optional) logical, true to suppress the
%           warnings of the initialization;
%
% OUTPUT
% - PDE_out:    the initialized PDE: a 'pde_struct' object for the terms
%           format, a 'struct' for the batch and legacy terms formats; an
%           empty 'pde_struct' object if no input is given.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ, 08/18/2022.
% MMP, 11/19/2026: Release comment pass (comment-only).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==0
    PDE_out = pde_struct();
    return
end
    
% Call the appropriate initialization function.
if isa(PDE,'pde_struct') || isa(PDE,'struct') && isfield(PDE,'x')
    % 'pde_struct', or 'struct' with field x: pde_struct initialization
    if isa(PDE,'struct')
        PDE = pde_struct(PDE);
    end
    if nargin==1
        PDE_out = initialize(PDE);
    else
        PDE_out = initialize(PDE,silent_initialize);
    end
% Batch format (fields n0, n1, n2)
elseif isfield(PDE,'n0') || isfield(PDE,'n1') || isfield(PDE,'n2')
    PDE_out = initialize_PIETOOLS_PDE_batch(PDE);
% Terms format (fields n, ODE or PDE); silent through the base-workspace
% flag silent_initialize_pde
elseif isfield(PDE,'n') || isfield(PDE,'ODE')  || isfield(PDE,'PDE')
    if nargin==2 && silent_initialize
        evalin('base','silent_initialize_pde=true;');
    end
    PDE_out = initialize_PIETOOLS_PDE_terms_legacy(PDE);
else
    error('The input PDE is not appropriately specified. Please define your PDE as a "pde_struct" class object, and consult the manual and examples for illustration of the structure.')
end

end