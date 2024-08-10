function [PDE_out] = vertcat(varargin)
% PDE_OUT = VERTCAT(PDE_1,PDE_2) vertically concatenates the PDEs defined
% by PDE_1 and PDE_2.
% 
% INPUT
% - PDE_1, PDE_2, PDE_3, ...:   
%                   'pde_struct' objects, either all representing PDE
%                   equations (for different state variables and outputs),
%                   or all representing free terms to be used to
%                   construct PDEs.
%
% OUTPUT
% - PDE_out:        'pde_struct' object containing the equations from all
%                   input PDE structures, concatenated such that the
%                   equations from PDE_1 appear first, those from PDE_2
%                   appear second, etc..
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
% % If more than two inputs are provided, just vertcat repeatedly.
if nargin==1
    % Vertcat of a single object is just the object;
    PDE_out = varargin{1};
    return
elseif nargin==2
    % Vertcat of two objects we process below;
    PDE_1 = varargin{1};    PDE_2 = varargin{2};
else
    % Vertcat of more than two objects we just keep repeating;
    PDE_1 = varargin{1};    PDE_2 = varargin{2};
end

% % Check that inputs are of appropriate type.
% We only support concatenation of PDE structures with other PDE
% structures, or with zeros, nothing else.
if isa(PDE_1,'state')
    PDE_1 = state2pde_struct(PDE_1);
elseif isa(PDE_1,'terms')
    PDE_1 = terms2pde_struct(PDE_1);
end
if isa(PDE_1,'polynomial') && isdouble(PDE_1)
    % Just in case, for some reason, zeros are specified as polynomial...
    PDE_1 = double(PDE_1);
end
if ~isa(PDE_1,'pde_struct') && (~isnumeric(PDE_1) || ~all(all(PDE_1==0)))
    error("Concatenation of 'pde_struct' objects with non-'pde_struct' objects is not supported.")
end
if isa(PDE_2,'state')
    PDE_2 = state2pde_struct(PDE_2);
elseif isa(PDE_2,'terms')
    PDE_2 = terms2pde_struct(PDE_2);
end
if isa(PDE_2,'polynomial') && isdouble(PDE_2)
    PDE_2 = double(PDE_2);
end
if ~isa(PDE_2,'pde_struct') && (~isnumeric(PDE_2) || ~all(all(PDE_2==0)))
    error("Concatenation of 'pde_struct' objects with non-'pde_struct' objects is not supported.")
end

% % Convert single PDE variables (states, inputs, outputs) to PDE terms.
is_pde_term_1 = false;      % Keep track of whether PDE_2 is a free term.
if ~isnumeric(PDE_1)
    [is_pde_var_1,obj] = is_pde_var(PDE_1);
    if is_pde_var_1
        PDE_1 = var2term(PDE_1,obj);
        is_pde_term_1 = true;
    else
        is_pde_term_1 = is_pde_term(PDE_1);
    end
else
    is_pde_term_1 = true;   % Treat 0 as a term.
end
is_pde_term_2 = false;      % Keep track of whether PDE_2 is a free term.
if ~isnumeric(PDE_2)
    [is_pde_var_2,obj] = is_pde_var(PDE_2);
    if is_pde_var_2
        PDE_2 = var2term(PDE_2,obj);
        is_pde_term_2 = true;
    else
        is_pde_term_2 = is_pde_term(PDE_2);
    end
else
    is_pde_term_2 = true;   % Treat 0 as a term.
end

% % Make sure either both inputs correspond to completed PDEs,
% % or both inputs correspond to loose terms to add to PDEs
if (~is_pde_term_1 && is_pde_term_2)
    % Convert free terms to equation.
    PDE_2 = (PDE_2==0);
elseif (is_pde_term_1 && ~is_pde_term_2)
    % Convert free terms to equation.
    PDE_1 = (PDE_1==0);
    %error("Concatenation of completed PDE equations with separate PDE terms or zeros is not supported.")
end


% % % Perform the actual concatenation
% % First, combine the state components, inputs, and outputs of the two
% % PDEs, so that both depend are expressed in terms of the same objects.
if ~isnumeric(PDE_1) && ~isnumeric(PDE_2)
    [PDE_1,PDE_2] = pde_common_basis(PDE_1,PDE_2);
end

% % Consider four cases:
if isnumeric(PDE_1)
    % % Case 1: Concatenate zeros with PDE_2;
    PDE_out = PDE_2;
    PDE_out.free = [cell(size(PDE_1,1),1); PDE_2.free];
elseif isnumeric(PDE_2)
    % % Case 1: Concatenate PDE_1 with zeros;
    PDE_out = PDE_1;
    PDE_out.free = [PDE_1.free; cell(size(PDE_2,2))];
elseif is_pde_term(PDE_1)
    % % Case 3: Both PDE structures correspond to loose terms, yet to be
    % % used to define an actual PDE. In this case, just concatenate
    % % elements of PDE.free, defining these loose terms.
    PDE_out = PDE_1;
    PDE_out.free = [PDE_1.free; PDE_2.free];
else
    % % Case 4: Both PDE structures correspond to alreay completed PDE
    % % and output equations. In this case, concatenate equations of each
    % % type: 'x', 'y', 'z', and 'BC'.
    % % Note that "pde_common_basis" already concatenated the fields 'x',
    % % 'y', and 'z', so all that needs to be done there is just add the
    % % terms defining the associated equation.
    PDE_out = PDE_1;
    eq_types = {'x';'y';'z'};
    for kk=1:numel(eq_types)
        % % Since PDE_out is initialized as PDE_1, make sure to add
        % % equations from PDE_2, where appropriate.
        obj = eq_types{kk};
        for ii=1:numel(PDE_out.(obj))
            if ~isfield(PDE_2.(obj){ii},'term') || isempty(PDE_2.(obj){ii}.term)
                % PDE_2 does not specify an equation for this object
                % --> stick with the equation from PDE_1.
                continue
            elseif ~isfield(PDE_1.(obj){ii},'term') || isempty(PDE_1.(obj){ii}.term)
                % PDE 1 does not specify an equation for this object, but
                % PDE 2 does
                % --> copy equation from PDE_2.
                PDE_out.(obj){ii}.term = PDE_2.(obj){ii}.term;
                if isfield(PDE_2.(obj){ii},'tdiff')
                    PDE_out.(obj){ii}.tdiff = PDE_2.(obj){ii}.tdiff;
                end
            else
                % Both PDE structures specify an equation for this object,
                % this may lead to conflict...
                error(["The two PDEs to concatenate cannot both specify an equation for the same object '",obj,"'."])
            end
        end
    end
    % Concatenate the boundary conditions.
    PDE_out.BC = [PDE_1.BC; PDE_2.BC];
    % Note that the returned PDE is no longer initialized.
    PDE_out.is_initialized = false;
end

% % Deal with remaining PDEs to concatenate.
if nargin>=3
    PDE_out = vertcat(PDE_out,varargin{3:end});
end

end