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
% DJ, 01/06/2025: Expand support for equating vector-valued terms;
% DJ, 08/22/2025: Allow equations x==0 or 0==x to be set, adjusting size of
%                   0 to match that of x;

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
if ~isnumeric(LHS)                                                          % DJ, 01/06/2025
    % Make sure that the vector-valued size of each row of terms in LHS
    % matches that in RHS
    sz_L = zeros(numel(LHS.free),1);
    for ii=1:numel(LHS.free)
        if isempty(LHS.free{ii})                                            % DJ, 08/22/2025
            % Assume LHS.free{ii} = zeros(sz_R(ii));
            sz_L(ii) = nan;
        else
            sz_L(ii) = LHS.free{ii}.size;
        end
    end
    sz_R = zeros(numel(RHS.free),1);
    for ii=1:numel(RHS.free)
        if isempty(RHS.free{ii})                                            % DJ, 08/22/2025
            % Assume RHS.free{ii} = zeros(sz_L(ii));
            sz_R(ii) = nan;
        else
            sz_R(ii) = RHS.free{ii}.size;
        end
    end
    if any(isnan(sz_L)) || any(isnan(sz_R))                                 % DJ, 08/22/2025
        % If there is only one uncertain size, match the LHS and RHS
        % dimensions
        if sum(isnan(sz_L))==1 && ~any(isnan(sz_R)) && (sum(sz_R)-sum(sz_L(~isnan(sz_L))))>=0
            sz_L(isnan(sz_L)) = sum(sz_R) - sum(sz_L(~isnan(sz_L)));
        elseif sum(isnan(sz_R))==1 && ~any(isnan(sz_L)) && (sum(sz_L) - sum(sz_R(~isnan(sz_R))))>=0
            sz_R(isnan(sz_R)) = sum(sz_L) - sum(sz_R(~isnan(sz_R)));
        else
            % Otherwise, if the number of equations match, get size of LHS 
            % from that of RHS, and vice versa
            if numel(sz_L)==numel(sz_R) && ~any(isnan(sz_L) & isnan(sz_R))
                sz_L(isnan(sz_L)) = sz_R(isnan(sz_L));
                sz_R(isnan(sz_R)) = sz_L(isnan(sz_R));
            else
                error("Equality cannot be set: size of zeros on LHS or RHS is ambiguous.")
            end
        end
    end
    if sum(sz_L)~=sum(sz_R)
        error("Number of rows of terms to equate does not match; equality cannot be enforced.")
    elseif numel(sz_L)~=numel(sz_R) || ~all(sz_L==sz_R)
        % It seems we cannot just set each row of terms in LHS equal to a
        % corresponding row of terms in RHS. Instead, we will combine all
        % rows of terms in both LHS and RHS into a single vector-valued row
        % of terms, and enforce only a single vector-value equality.
        LHS = squash_eqs(LHS,sz_L);
        RHS = squash_eqs(RHS,sz_R);
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
    if isempty(RHS.free{ii})                                                % DJ, 08/22/2025
        % Remove equations 0==0
        continue
    end
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
        if isfield(RHS.free{ii}.term{1},'C') && ~all(all(isequal(RHS.free{ii}.term{1}.C,eye(PDE_out.(obj){eq_num}.size))))
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



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function PDE_new = squash_eqs(PDE,sz_eqs)
% PDE_NEW = SQUASH_EQS(PDE) takes a PDE structure describing free terms (no
% completed equations), and squashes all rows of terms (each row to be
% added to a separate equation) into a single vector-valued row of terms.
%
% INPUT
% - PDE:        'pde_struct' object representing a set of n rows of free
%               terms to be used to define equations, specified through the
%               fields PDE.free{ii}.term, for ii=1:n;
% - sz_eqs:     nx1 array specifying for each row of terms the vector size
%               of those terms, sz_eqs(ii) = PDE.free{ii}.size;
%
% OUTPUT
% - PDE_new:    'pde_struct' object representing the same equations, but
%               now with PDE_new.free only containing a single element,
%               with PDE_new.free{1}.term including all terms from all
%               elements PDE.free{ii}, premultiplied with a suitable
%               permutation matrix;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 01/06/2025: Initial coding;
%

% Initialize a new object with just 1 vector-valued equation
PDE_new = PDE;
PDE_new.free = PDE.free(1);
% Augment all terms in first row to new size by premultiplying coefficients.
C_aug = [eye(sz_eqs(1));zeros(sum(sz_eqs(2:end)),sz_eqs(1))];
for jj=1:numel(PDE.free{1}.term)
    if isfield(PDE_new.free{1}.term{jj},'C')
        PDE_new.free{1}.term{jj}.C = C_aug*PDE_new.free{1}.term{jj}.C;
    else
        PDE_new.free{1}.term{jj}.C = C_aug;
    end
end
% Keep track of the variables that appear in the first row;
var_list = cell(size(PDE.free{1}.vars,1),1);
for kk=1:numel(var_list)
    var_list{kk} = PDE.free{1}.vars(kk).varname{1};
end
% Add all other terms, augmented to new size.
for ii=2:numel(sz_eqs)
    trm_cntr = numel(PDE_new.free{1}.term);
    PDE_new.free{1}.term = [PDE_new.free{1}.term,PDE.free{ii}.term];
    C_aug = [zeros(sum(sz_eqs(1:ii-1)),sz_eqs(ii));
             eye(sz_eqs(ii));
             zeros(sum(sz_eqs(ii+1:end)),sz_eqs(ii))];
    for jj=trm_cntr+(1:numel(PDE.free{ii}.term))
        if isfield(PDE_new.free{1}.term{jj},'C')
            PDE_new.free{1}.term{jj}.C = C_aug*PDE_new.free{1}.term{jj}.C;
        else
            PDE_new.free{1}.term{jj}.C = C_aug;
        end
    end
    % Keep track of the variables that appear in all rows.
    for kk=1:size(PDE.free{ii}.vars,1)
        varname_kk = PDE.free{ii}.vars(kk).varname{1};
        if ~ismember(varname_kk,var_list)
            var_list = [var_list;{varname_kk}];
        end
    end
end
% Set the variables and size of the new squashed row of terms
PDE_new.free{1}.size = sum(sz_eqs);
PDE_new.free{1}.vars = polynomial(var_list);
% Continue with new RHS
PDE = PDE_new;

end