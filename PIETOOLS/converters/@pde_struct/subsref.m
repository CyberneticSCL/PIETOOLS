function val = subsref(PDE,prop)
% SUBSREF allows properties and elements of the "pde_struct" object PDE to
% be retrieved. The function is automatically called when calling
% e.g. PDE(i,j) or PDE.obj for a "pde_struct" PDE.
%
% INPUTS:
% - PDE:    A pde_struct object of which to extract a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PDE
%           structure to retrieve.
%
% OUTPUTS:
% - val:    The value of the field of "PDE" specified by "prop".
%
% EXAMPLE: calling PDE.x{1}.term{2}.C(3), we have
%   PDE_in = PDE;
%   prop(1).type = '.',     prop(1).subs = 'x';
%   prop(2).type = '{}',    prop(2).subs = {[1]};
%   prop(3).type = '.',     prop(3).subs = 'term';
%   prop(4).type = '{}',    prop(4).subs = {[2]};
%   prop(5).type = '.',     prop(5).subs = 'C';
%   prop(6).type = '()',    prop(6).subs = {[3]};
%
% NOTES:
% For the most part, this function just calls the Matlab built-in subsref
% function. However, in the display, certain coefficients are shown as
% C_{ij}, so we added a feature to extract these coefficients by calling
% PDE.C{i,j}.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022 PIETOOLS Team
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
% DJ, 11/29/2022: Initial coding;
% DJ, 01/03/2025: Update to assume a loose PDE variable is specified as a
%                   single free term, see also update to "pde_var";
% DJ, 01/04/2025: Allow order of differentiability of state to be extracted;


% % % If the object is a single PDE variable, we have a separate case for
% % % extracting information.
[is_var,obj] = is_pde_var(PDE);

% At this point, we just use the built-in subsasgn function. Additional
% features such as '()' type subsasgn may be introduced in later updates.
if strcmp(prop(1).type,'{}')
    error("''{}'' type subsref is currently not supported for 'pde_struct' class objects.")
elseif strcmp(prop(1).type,'()')
    % % Allow terms to be retrieved by calling PDE(i,j);
    % % Process the inputs
    sz_PDE = size(PDE);
    if isscalar(prop(1).subs) && sz_PDE(1)==1
        % For single input case, if there is only one equation, 
        % extract the desired term.
        prop(1).subs = {1,prop(1).subs{1}};
    elseif isscalar(prop(1).subs)
        % In single input case for multiple equations, 
        % assume full equations are desired.
        prop(1).subs = {prop(1).subs{1},':'};
    elseif numel(prop(1).subs)>2
        error("Subsref is only supported with two inputs: equation number and term number.")
    end
    % Extract equation and term numbers
    eq_idcs = prop(1).subs{1};       term_idcs = prop(1).subs{2};
    if ischar(eq_idcs) && strcmp(eq_idcs,':')
        % Assume all equations are desired.
        eq_idcs = 1:sz_PDE(1);
    elseif ~isnumeric(eq_idcs) || length(eq_idcs)~=numel(eq_idcs)
        error("Equation indices should be specified as 1xn array of nonnegative integers.")
    elseif isnumeric(eq_idcs)
        if max(eq_idcs)>sz_PDE(1)
            error("Equation indices cannot exceed the number of equations in the PDE!")
        elseif any(eq_idcs<=0) || any(round(eq_idcs)~=eq_idcs)
            error("Equation indices should be specified as 1xn array of nonnegative integers.")
        end
    end
    if ischar(term_idcs) && strcmp(term_idcs,':')
        % Assume all equations are desired.
        term_idcs = 1:sz_PDE(2);
    elseif ~isnumeric(term_idcs) || length(term_idcs)~=numel(term_idcs)
        error("Equation indices should be specified as 1xn array of integers.")
    elseif isnumeric(eq_idcs)
        if max(term_idcs)>sz_PDE(2)
            error("Term indices cannot exceed the maximal number of terms that appear in the equations!")
        end
    end

    % % Build the output
    val = PDE;
    if ~isempty(PDE.free)
        % % The input PDE describes only free terms
        % % --> just return the desired terms
        free_tmp = PDE.free;
        for ii=1:numel(free_tmp)
            if isfield(free_tmp{ii},'term') && ~isempty(free_tmp{ii}.term)
                % Retain only terms specified by term_idcs.
                term_idcs_ii = term_idcs(term_idcs<=numel(free_tmp{ii}.term));
                free_tmp{ii}.term = free_tmp{ii}.term(term_idcs_ii);
            end
        end
        val.free = free_tmp(eq_idcs);
    else
        % % The input PDE describes actual equations
        % % --> we'll have to move the terms into the field 'free'.
        % Equations can be extracted of each of the following types.
        objs_arr = {'x';'y';'z';'BC'};
        % Get the number of equations of each type.
        neqs_arr = size(PDE,objs_arr);
        nneqs_arr = cumsum([0;neqs_arr]);
        % Initiale a cell of free equations for the output.
        free_tmp = cell(length(eq_idcs),1);
        % Allow the left-hand side of the equation to be extracted as well.
        if any(term_idcs==0)
            if numel(term_idcs)~=1
                error("When extracting term number 0, nu other terms can be extracted.")
            % elseif any(eq_idcs)>nneqs_arr(end-1)
            %     error("Left-hand side of boundary conditions cannot be extracted.")
            end
        end
        for ii=1:numel(free_tmp)
            % Determine which object equation ii corresponds to
            obj_idx = find(eq_idcs(ii)<=nneqs_arr,1,'first')-1;
            obj = objs_arr{obj_idx};
            eq_num = eq_idcs(ii)-nneqs_arr(obj_idx);
            % Extract size and variable information on the equation
            eq_struct = PDE.(obj){eq_num};
            if isfield(eq_struct,'size')
                free_tmp{ii}.size = eq_struct.size;
            end
            if isfield(eq_struct,'vars')
                free_tmp{ii}.vars = eq_struct.vars;
            end
            free_tmp{ii}.term = {};
            % Extract the desired terms from the equation
            if all(term_idcs==0)
                % Extract left-hand side of equation.
                if obj_idx==4
                    % Left-hand side of boundary condition is zero.
                    continue
                else
                    % Left-hand side of output or PDE is just the output or
                    % state.
                    free_tmp{ii}.term{1}.(obj) = eq_num;
                    free_tmp{ii}.term{1}.C = eye(PDE.(obj){eq_num}.size);
                    if isfield(eq_struct,'tdiff')
                        free_tmp{ii}.term{1}.tdiff = eq_struct.tdiff;
                    end
                end
            else
                if isfield(eq_struct,'term') && ~isempty(eq_struct.term)
                    % Retain only terms specified by term_idcs.
                    term_idcs_ii = term_idcs(term_idcs<=numel(eq_struct.term));
                    free_tmp{ii}.term = eq_struct.term(term_idcs_ii);
                end
            end
        end
        % % Remove all non-free equations.
        for ii=1:nneqs_arr(end-1)
            % Determine which object equation ii corresponds to
            obj_idx = find(ii<=nneqs_arr,1,'first')-1;
            obj = objs_arr{obj_idx};
            eq_num = ii-nneqs_arr(obj_idx);
            % Remove all terms from the equation.
            % if isfield(val.(obj){eq_num},'term')
            %     val.(obj){eq_num} = rmfield(val.(obj){eq_num},'term');
            % end%
            val.(obj){eq_num}.term = {};
        end
        val.BC = cell(0,1);
        % % Finally, set the free equations in the PDE.
        val.free = free_tmp;
        val.is_initialized = false;
    end
    % Perform remaining subsref.
    if numel(prop)>1
        val = subsref(val,prop(2:end));
    end

elseif strcmp(prop(1).type,'.')
    % % % Allow,size, variables and domain of PDE variable to be extracted.
    if is_var && (strcmp(prop(1).subs,'size') || strcmp(prop(1).subs,'vars') || strcmp(prop(1).subs,'var') || strcmp(prop(1).subs,'dom') || strcmp(prop(1).subs,'diff'))
        if strcmp(prop(1).subs,'var')
            prop(1).subs = 'vars';
        end
        if strcmp(prop(1).subs,'diff') && ~isfield(PDE.(obj){1},'diff')
            val = zeros(1,0);
        else
            val = PDE.(obj){1}.(prop(1).subs);
        end
        
        % Perform remaining subsref.
        if numel(prop)>1
            val = builtin('subsref',val,prop(2:end));
        end
        return
    end
    % % % Allow coefficients to be retrieved by calling PDE.C{i,j}
    if strcmpi(prop(1).subs,'C')
        % Deal with case that the PDE corresponds to a single term.
        sz_PDE = size(PDE);
        if all(sz_PDE==1) && (isscalar(prop) || strcmp(prop(2).type,'()'))
            % Allow coefficient C to be extracted in case of a single term.
            prop_new = prop(1);
            prop_new(2).type = '{}';        % Extract coefficients {1,1};
            prop_new(2).subs = {1,1};
            prop_new = [prop_new,prop(2:end)];
            prop = prop_new;
        elseif isscalar(prop)
            error('To extract the coefficients C_{ij} appearing in term j of equation i, call "PDE.C{i,j}".')
        end
        if ~strcmp(prop(2).type,'{}')
            error('To extract the coefficients C_{ij} appearing in term j of equation i, call "PDE.C{i,j}".')
        elseif numel(prop(2).subs)<2
            error('Both an equation number i and term number j must be specified to extract coefficients "PDE.C{i,j}".')
        elseif numel(prop(2).subs)>3
            error('Only an equation number i, term number j, and factor number k can be specified when extracting coefficients "PDE.C{i,j,k}".')
        else
            ii = prop(2).subs{1};
            jj = prop(2).subs{2};
            if numel(prop(2).subs)==3
                kk = prop(2).subs{3};
            else
                kk = 1;
            end
        end
        % Check that the equation and term number make sense.
        if ii <= 0 || ii~=round(ii) || ii~=real(ii)
            error('Equation number i must be specified as real, nonnegative integer when extracting coefficients "PDE.C{i,j}".')
        elseif jj <= 0 || jj~=round(jj) || jj~=real(jj)
            error('Term number j must be specified as real, nonnegative integer when extracting coefficients "PDE.C{i,j}".')
        % elseif ii>sz_PDE(1)
        %     error("Equation number i cannot exceed number of equations in the PDE, when extracting coefficients 'PDE.C{i,j}'.")
        elseif jj>sz_PDE(2)
            error("Term number j cannot exceed maximum number of terms in each equation, when extracting coefficients 'PDE.C{i,j}'.")
        elseif kk <= 0 || kk~=round(kk) || kk~=real(kk)
            error('Term number j must be specified as real, nonnegative integer when extracting coefficients "PDE.C{i,j}".')
        end
        % Distinguish between free equation, dynamics, outputs, or BCs
        if ~isempty(PDE.free)
            eq = PDE.free{ii};
            eqname = ['free{',num2str(ii),'}'];
        elseif ii <= numel(PDE.x)
            % Coefficients in PDE dynamics are requested.
            eq = PDE.x{ii};
            eqname = ['x{',num2str(ii),'}'];
        elseif ii <= numel(PDE.x) + numel(PDE.y)
            % Coefficients in observed output equations are requested.
            eqnum = ii - numel(PDE.x);
            eq = PDE.y{eqnum};
            eqname = ['y{',num2str(eqnum),'}'];
        elseif ii <= numel(PDE.x) + numel(PDE.y) + numel(PDE.z)
            % Coefficients in regulated output equations are requested.
            eqnum = ii - numel(PDE.x) - numel(PDE.y);
            eq = PDE.z{eqnum};
            eqname = ['z{',num2str(eqnum),'}'];
        elseif ii <= numel(PDE.x) + numel(PDE.y) + numel(PDE.z) + numel(PDE.BC)
            % Coefficients in BCs are requested.
            eqnum = ii - numel(PDE.x) - numel(PDE.y) - numel(PDE.z);
            eq = PDE.BC{eqnum};
            eqname = ['BC{',num2str(eqnum),'}'];
        else
            error('The proposed equation index exceeds the number of equations in the system.')
        end
        % Extract coefficients from the desired term.
        if jj <= 0 || jj~=round(jj) || jj~=real(jj)
            error('Term number j must be specified as real, nonnegative integer when extracting coefficients "PDE.C{i,j}".')
        elseif jj>numel(eq.term)
            error(['The proposed term index exceeds the number of terms in the equation "PDE.',eqname,'".'])
        else
            trm = eq.term{jj};
        end
        if kk>numel(trm)
            error("Factor number k cannot exceed the number of factors in the desired term.")
        else
            val = trm(kk).C;
        end
        % Perform remaining subsref.
        if numel(prop)>2
            val = builtin('subsref',val,prop(3:end));
        end
    else
        % Otherwise, we have nothing fancy implemented at this time.
        val = builtin('subsref',PDE,prop);
    end
end

end