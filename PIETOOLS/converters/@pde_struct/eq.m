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
% Copyright (C)2024 PIETOOLS Team
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
% DJ, 12/28/2024: Add support for concatenating single vector LHS with
%                 RHS consisting of multiple scalar-valued equations; 
% DJ, 12/29/2024: Bugfix minus sign in case LHS has coefficients.
% DJ, 01/03/2025: Use "is_zero" to indicate zero equations (e.g. y==0);
% DJ, 01/05/2025: Bugfix, replace {eq_num} with {eq_num,1};

% % % Process the inputs

% % Check that the first input makes sense
if isa(LHS,'state')
    LHS = state2pde_struct(LHS);
elseif isa(LHS,'terms')
    LHS = terms2pde_struct(LHS);
end
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
    error("For declaring an equation, the left-hand side should correspond to the temporal derivative of a state, or to an output signal.")
end

% % Check that the second input makes sense
if isa(RHS,'state')
    RHS = state2pde_struct(RHS);
elseif isa(RHS,'terms')
    RHS = terms2pde_struct(RHS);
end
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
    error("Right-hand side in equation of 'pde_struct' objects should consist solely of PDE terms.")
end

% % Check that the number of rows of LHS and RHS match
if ~isnumeric(LHS) && numel(LHS.free)~=numel(RHS.free)
    % The number of objects does not match, but maybe their sizes add
    % up to match. Currently, we only support this if LHS describes a 
    % single (vector-valued) object.
    if isscalar(LHS.free)
        sz_L = LHS.free{1}.size;
        nR = numel(RHS.free);
        sz_R = RHS.free{1}.size;
        vars_R = RHS.free{1}.vars;

        % Make sure the sizes and variables of all objects match.
        for ii=2:nR
            if isfield(RHS.free{ii},'size') && RHS.free{ii}.size~=sz_R
                error("For equating vector-valued object with multiple terms, all terms must have same size.")
            elseif isfield(RHS.free{ii},'vars')
                if ~(isempty(RHS.free{ii}.vars) && isempty(vars_R)) && ...
                    (~all(size(RHS.free{ii}.vars)==size(vars_R)) || ~all(isequal(RHS.free{ii}.vars,vars_R)))
                    error("For equating vector-valued object with multiple terms, all terms must have same variables.")
                end
            end
        end
        if nR*sz_R~=sz_L
            error("Number of objects to equate does not match; concatenation not supported.")
        end
        
        % Initialize a new object with just 1 vector-valued
        RHS_new = RHS;
        RHS_new.free = RHS.free(1);
        RHS_new.free{1}.size = sz_L;
        % Augment all terms in RHS.free{1} to new size.
        C_aug = [eye(sz_R);zeros((nR-1)*sz_R,sz_R)];
        for jj=1:numel(RHS.free{1}.term)
            if isfield(RHS_new.free{1}.term{jj},'C')
                RHS_new.free{1}.term{jj}.C = C_aug*RHS_new.free{1}.term{jj}.C;
            else
                RHS_new.free{1}.term{jj}.C = C_aug;
            end
        end
        % Add all other terms, augmented to new size.
        for ii=2:nR
            cntr = numel(RHS_new.free{1}.term);
            RHS_new.free{1}.term = [RHS_new.free{1}.term,RHS.free{ii}.term];
            C_aug = [zeros((ii-1)*sz_R,sz_R);eye(sz_R);zeros((nR-ii)*sz_R,sz_R)];
            for jj=cntr+(1:numel(RHS.free{ii}.term))
                if isfield(RHS_new.free{1}.term{jj},'C')
                    RHS_new.free{1}.term{jj}.C = C_aug*RHS_new.free{1}.term{jj}.C;
                else
                    RHS_new.free{1}.term{jj}.C = C_aug;
                end
            end
        end
        % Continue with new RHS
        RHS = RHS_new;
    else
        error("Number of objects to equate does not match; concatenation currently not supported.")
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
            PDE_out.(obj){eq_num,1}.term = RHS.free{ii}.term(2:end);
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
            PDE_out.(obj){eq_num,1}.tdiff = tdiff;
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
                    PDE_out.(obj){eq_num,1}.term{jj}.C = -inv(C_LHS)*PDE_out.(obj){eq_num}.term{jj}.C;    % DJ, 12/29/2024
                else
                    PDE_out.(obj){eq_num,1}.term{jj}.C = -1;
                end
            end
        else
            % Take negative coefficients to account for fact that terms
            % have been moved
            for jj=1:numel(PDE_out.(obj){eq_num}.term)
                if isfield(PDE_out.(obj){eq_num}.term{jj},'C')
                    PDE_out.(obj){eq_num,1}.term{jj}.C = -PDE_out.(obj){eq_num}.term{jj}.C;
                else
                    PDE_out.(obj){eq_num,1}.term{jj}.C = -1;
                end
            end
        end
        % If there are no other terms, we must have e.g. d/dt x = 0.
        if isscalar(RHS.free{ii}.term)                                      % DJ, 01/03/2025
            PDE_out.(obj){eq_num,1}.is_zero = true;
        end
    else
        % % The equation corresponds to a boundary condition.
        BC_num = numel(PDE_out.BC)+1;
        PDE_out.BC{BC_num,1}.term = RHS.free{ii}.term(1:end);
    end
    % Get rid of loose PDE terms.
    PDE_out.free = {};
end

end