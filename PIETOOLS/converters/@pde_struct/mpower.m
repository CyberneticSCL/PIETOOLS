function PDE_out = mpower(PDE_in,deg)
% PDE_OUT = POWER(PDE_IN,DEG) computes the power of a variable or term in 
% a PDE.
%
% INPUTS
% - PDE_in: 'pde_struct' object representing a single term or set of terms 
%           in a PDE. Can be at most 1 row of scalar-valued terms;
% - deg:    scalar integer specifying to which power to raise the terms in
%           the PDE;
%
% OUTPUTS
% - PDE_out:    'pde_struct' object representing the power of the terms in
%               the input PDE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 02/18/2026: Initial coding;

% Make sure the inputs are of appropriate type
if ~isa(PDE_in,'pde_struct')
    error("The first inputs must be a 'pde_struct' object.")
end
if ~isa(deg,'double') || ~isscalar(deg) || deg<0 || deg~=round(deg)
    error("Degrees must be specified as real nonnegative integer.")
end

% Make sure the input corresponds to a single state/input/output
% or set of terms
[is_pde_var_in,obj_in] = is_pde_var(PDE_in);
if (~is_pde_var_in && ~is_pde_term(PDE_in))
    error("Equations cannot be raised to a power.")
end
% Convert the PDE variable to a term
if is_pde_var_in
    PDE_in = var2term(PDE_in,obj_in);
end

% Make sure the structure defines a single row of scalar-valued terms
nr = size(PDE_in,'free','vec_size_tot');
if nr>1
    error("Taking power of vector-valued variables is not supported. Use '.^' instead for computing the elementwise power.")
end

PDE_out = power(PDE_in,deg);

end