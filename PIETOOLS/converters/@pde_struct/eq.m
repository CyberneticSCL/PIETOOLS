function PDE_out = eq(LHS,RHS)
% PDE_OUT = EQ(LHS,RHS) declares a PDE equation setting LHS equal to RHS.
%
% INPUT
% - LHS:    pde_struct object representing the left-hand side of the
%           equation, either a temporal derivative of a state x, or an
%           output y or z. Can also be 0 to declare boundary conditions.
% - RHS:    pde_struct object representing the right-hand side of the
%           equation, corresponding to terms in the PDE or output equation.
%           Can also be set to 0 to declare boundary conditions.
%
% OUTPUT
% - PDE_out:    pde+struct object representing the PDE or output equation
%               LHS=RHS;
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

% % Check that the first input makes sense
if isa(LHS,'polynomial') && isdouble(LHS)
    LHS = double(LHS);
elseif isa(LHS,'polynomial')
    error("Equating PDE objects with polynomials is not supported.")
end
if isnumeric(LHS) && ~all(all(LHS==0))
    error("Equating PDE objects with nonzero values is not supported.")
elseif ~isnumeric(LHS) && ~isa(LHS,'pde_struct')
    error("Equating 'pde_struct' objects with non-'pde_struct' objects is not supported.")
elseif isa(LHS,'pde_struct') && ~is_pde_term(LHS)
    if is_pde_var(LHS)
        LHS = var2term(LHS);
    else
        error("For declaring an equation, the left-hand side should correspond to the temporal derivative of a state, or to an output signal.")
    end
end

% % Check that the second input makes sense
if isa(RHS,'polynomial') && isdouble(RHS)
    RHS = double(RHS);
elseif isa(RHS,'polynomial')
    error("Equating PDE objects with polynomials is not supported.")
end
if isnumeric(RHS) && ~all(all(RHS==0))
    error("Equating PDE objects with nonzero values is not supported.")
elseif isnumeric(RHS)
    % Move zeros for BCs to left-hand side.
    RHS = LHS;
    LHS = 0;
elseif ~isa(RHS,'pde_struct')
    error("Equating 'pde_struct' objects with non-'pde_struct' objects is not supported.")
elseif ~is_pde_term(RHS)
    if is_pde_var(RHS)
        RHS = var2term(RHS);
    else
        error("Right-hand side in equation of 'pde_struct' objects should consist solely of PDE terms.")
    end
end

% % Check that the number of rows of LHS and RHS match
if ~isnumeric(LHS)
    if numel(LHS.free)~=numel(RHS.free)
        error("Number of rows of terms to equate does not match.")
    end
end


% % % Proceed only with equations of the form 0==...
if ~isnumeric(LHS)
    RHS = minus(LHS,RHS);
    LHS = 0;
end

% % % Loop over all equations
PDE_out = RHS;
for ii=1:numel(RHS.free)
    % % Check if the first term corresponds to a "left-hand side object",
    % % so the temporal derivative of a state component, or an output.
    [is_LHS_ii,obj,eq_num,tdiff] = is_LHS_term(RHS.free{ii}.term{1});
    if is_LHS_ii
        % % The equation corresponds to a PDE (d/dt x=...) or an output
        % % equation (y=..., z=...);
        % % In particular, equation "eq_num" of type "obj".
        % Add the remaining terms to this equation.
        if ~isfield(PDE_out.(obj){eq_num},'term') || isempty(PDE_out.(obj){eq_num}.term)
            PDE_out.(obj){eq_num}.term = RHS.free{ii}.term(2:end);
        else
            error("Multiple equations are specified for the same state or output variable; this is not supported.")
        end
        
        % Check that the size of the terms makes sense for the object.
        if PDE_out.(obj){eq_num}.size~=RHS.free{ii}.size
            error("The number of rows in the PDE terms to equate do not match.")
        end
        % Check that the variables appearing in the terms make sense for
        % the object.
        if any(~ismember(RHS.free{ii}.vars.varname,PDE_out.(obj){eq_num}.vars.varname))
            error("The equation cannot be set: the left-hand side and right-hand side depend on different spatial variables.")
        end

        % For state variables, include the order of the temporal
        % derivaitve.
        if strcmp(obj,'x')
            PDE_out.(obj){eq_num}.tdiff = tdiff;
        end
        % Check if coefficients have been specified for the left-hand side
        if isfield(RHS.free{ii}.term{1},'C') && ~all(all(RHS.free{ii}.term{1}.C==eye(PDE_out.(obj){eq_num}.size)))
            C_LHS = RHS.free{ii}.term{1}.C;
            if isa(C_LHS,'polynomial') && ~isdouble(C_LHS)
                error("Equation cannot be set: an output signal or a temporal derivative of a state has been multiplied with a polynomial")
            else
                C_LHS = double(C_LHS);
            end
            % Divide coefficients in each term by those on the LHS.
            for jj=1:numel(PDE_out.(obj){eq_num}.term)
                if isfield(PDE_out.(obj){eq_num}.term{jj},'C')
                    PDE_out.(obj){eq_num}.term{jj}.C = inv(C_LHS)*PDE_out.(obj){eq_num}.term{jj}.C;
                else
                    PDE_out.(obj){eq_num}.term{jj}.C = -1;
                end
            end
        else
            % Take negative coefficients to account for fact that terms
            % have been moved
            for jj=1:numel(PDE_out.(obj){eq_num}.term)
                if isfield(PDE_out.(obj){eq_num}.term{jj},'C')
                    PDE_out.(obj){eq_num}.term{jj}.C = -PDE_out.(obj){eq_num}.term{jj}.C;
                else
                    PDE_out.(obj){eq_num}.term{jj}.C = -1;
                end
            end
        end
    else
        % % The equation corresponds to a boundary condition.
        BC_num = numel(PDE_out.BC)+1;
        PDE_out.BC{BC_num}.term = RHS.free{ii}.term(1:end);
    end
    % Get rid of loose PDE terms.
    PDE_out.free = {};
end

end