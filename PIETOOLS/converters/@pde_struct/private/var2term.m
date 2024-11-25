function PDE_trm = var2term(PDE_var,obj)
%VAR2TERM Convert a PDE variable to a free term in a PDE.
%
% INPUT
% - PDE_var:    pde_struct object representing a single state variable or
%               input;
% - obj:        (optional) char object specifying what type of object the
%               PDE variable is
%               ('x' for state, 'w' for disturbance, 'u' for controlled
%               input, 'y' for sensed output, 'z' for regulated output);
%
% OUTPUT
% - PDE_trm:    pde_struct object representing a term to be added to a PDE,
%               corresponding to just the given state component or input.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 06/23/2024

% % % Process the inputs
% Make sure the input is a pde_struct object.
if ~isa(PDE_var,'pde_struct')
    error("Function only takes 'pde_struct' objects as arguments.")
end
% Make sure the input corresponds to a single state/input/output.
if nargin<2
    % Check what type of object the PDE variable is.
    [is_pde_var_in,obj] = is_pde_var(PDE_var);
else
    % If type of variable is provided, assume it is certain that PDE_in is
    % indeed a PDE variable.
    is_pde_var_in = true;
end
if ~is_pde_var_in
    error("Input should correspond to a single state, input, or output variable.")
end

% % % Declare the term
% Initialize output object and empty equation.
PDE_trm = PDE_var;
vars = PDE_var.(obj){1}.vars;
PDE_trm.free{1}.size = PDE_var.(obj){1}.size;
PDE_trm.free{1}.vars = vars;

% Set the term, corresponding to just the object: no differentiation,
% integration, substitution or multiplication.
if strcmp(obj,'y') || strcmp(obj,'z')
    % Output signals get special treatment, as they can only appear on the
    % left-hand side of an equation.
    PDE_trm.free{1}.term{1}.(obj) = 1;
    PDE_trm.free{1}.term{1}.C = 1;
else
    PDE_trm.free{1}.term{1}.(obj) = 1;
    PDE_trm.free{1}.term{1}.C = eye(PDE_var.(obj){1}.size);
    if strcmp(obj,'x')
        PDE_trm.free{1}.term{1}.loc = vars';
        PDE_trm.free{1}.term{1}.D = zeros(1,size(vars,1));
    end
    PDE_trm.free{1}.term{1}.I = cell(size(vars,1),1);
end

end