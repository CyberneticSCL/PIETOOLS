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
        error("Right-hand side in eqaution of 'pde_struct' objects should consist solely of PDE terms.")
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




% 
% % % % Deal with case that LHS=0
% % % % --> declare BCs.
% if isnumeric(LHS)
%     PDE_out = RHS;
%     % Move each PDE temporary equation into a boundary condition.
%     for ii=1:numel(RHS.free)
%         PDE_out.BC{ii} = RHS.free{ii};
%     end
%     % Get rid of the temporary equations.
%     PDE_out.free = {};
%     return
% end
% 
% % % % Deal with case that LHS is [d/dt x; y; z];
% PDE_out = RHS;
% %ncomps_x = zeros(numel(LHS.x),1);
% 
% % Make sure the number of equations matches.
% if numel(LHS.free)~=numel(RHS.free)
%     error("The number of rows of the left- and right-hand sides should match.")
% end
% 
% % Make sure that at least one of sides contains a temporal derivative of a
% % state, or an output signal....
% 
% 
% % Determine what type of equation is being set.
% is_obj_list = [numel(LHS.x)>0; numel(LHS.y)>0; numel(LHS.z)>0];
% if ~any(is_obj_list)
%     error("The left-hand side should correspond to the temporal derivative of a state, or to an output signal.")
% elseif sum(is_obj_list)>1
%     error("The left-hand side can correspond to only one of the following: a state x, a regulated output z, or a sensed output y.")
% end
% obj_list = ['x';'y';'z'];
% obj = obj_list(is_obj_list);
% 
% % Make sure the number of equations in LHS and RHS match.
% if numel(LHS.(obj))~=numel(RHS.free)
%     error("The number of rows of the left- and right-hand sides should match.")
% end
% 
% % % Combine the object tables
% tab_1 = LHS.([obj,'_tab']);     tab_2 = RHS.([obj,'_tab']);
% if ~isempty(tab_2)
%     % % The object appears both in LHS and RHS (can really only be state x)
%     % % --> Combine information on the considered object from LHS and RHS
% 
%     % Create a unique table containing the index numbers and sizes of the
%     % object components that appear both in LHS and RHS
%     tab_full = [tab_1;tab_2];
%     [new_IDs,idcs1,idcs2] = unique(tab_full(:,1),'stable');   % new_IDs = tab_full(idcs1,1);      tab_full(:,1) = new_IDs(idcs2);
%     PDE_out.([obj,'_tab']) = tab_full(idcs1,:);
% 
%     % Initialize the associated object in the output structure.
%     n_obj_LHS = size(tab_1,1);      n_obj = length(new_IDs);
%     PDE_out.(obj) = cell(n_obj,1);
%     % Copy components primarily from LHS, as that has potential temporal
%     % derivaitve information.
%     is_LHS_comp = idcs1<=n_obj_LHS;     
%     PDE_out.(obj)(is_LHS_comp) = LHS.(obj)(idcs1(is_LHS_comp));
%     % Fill in the gaps corresponding to components that appear only in RHS.
%     PDE_out.(obj)(~is_LHS_comp) = RHS.(obj)(idcs1(~is_LHS_comp)-n_obj_LHS);
% 
%     % Keep track of how old component numbers relate to new ones
%     new_obnums_RHS = idcs2(n_obj_LHS+1:end);   % new_idx = new_obnums_RHS(old_idx);
% else
%     % % The object does not appear in RHS
%     % % --> just copy information from LHS.
%     PDE_out.([obj,'_tab']) = tab_1;
%     PDE_out.(obj) = LHS.(obj);
%     idcs1 = 1:numel(LHS.(obj));
% end
% 
% % % Loop over all elements on LHS, and set the associated equation.
% for ii=1:numel(LHS.(obj))
%     eq_num = idcs1(ii);
%     if strcmp(obj,'x') && (~isfield(LHS.(obj){eq_num},'tdiff') || LHS.(obj){eq_num}.tdiff==0)
%         error(["No temporal derivative is take of state component ",num2str(eq_num)," on the left-hand side of the equation."])
%     elseif isfield(LHS.(obj){eq_num},'term') && ~isempty(LHS.(obj){eq_num}.term)
%         error("The left-hand side contains terms, this is not supported.")
%     end
%     % Check that the size and variables match.
%     if LHS.(obj){eq_num}.size~=RHS.free{ii}.size
%         error("The size of the terms on the left-hand side should match that of the terms on the right-hand side.")
%     elseif any(~ismember(RHS.free{ii}.vars.varname,LHS.(obj){eq_num}.vars.varname))
%         error("The equation cannot be set: the left-hand side and right-hand side depend on different spatial variables.")
%     end
%     % % Add the terms to the equation
%     PDE_out.(obj){ii}.term = RHS.free{ii}.term;
%     if strcmp(obj,'x')
%     for jj=1:numel(PDE_out.(obj){ii}.term)
%         % Make sure the state number that appears in the term matches the
%         % order of the states in the output PDE structure.
%         if isfield(PDE_out.(obj){ii}.term{jj},'x')
%             PDE_out.(obj){ii}.term{jj}.x = new_obnums_RHS(PDE_out.(obj){ii}.term{jj}.x);
%         end
%     end
%     end
% end
% 
% % Get rid of the temporary equations.
% PDE_out.free = {};

end