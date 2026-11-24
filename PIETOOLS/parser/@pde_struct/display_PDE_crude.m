function [] = display_PDE_crude(PDE,name)
% DISPLAY_PDE_CRUDE(PDE,NAME) displays the input PDE in the Command Window,
% giving a symbolic overview of the specified equations/terms using
% UNICODE. This function is intended to display "unfinished"
% (un-initialized) pde_struct objects.
%
% INPUT
% - PDE:    'pde_struct' object, representing a PDE to be displayed.
%           Equations will be printed in the Command Window.
% - name:   (optional) char, the name of the variable, used in the hint
%           "name.C{i,j}" on extracting coefficients; 'PDE' by default.
%
% OUTPUT
% - none; the equations (or free terms) are printed, one per line.
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 06/23/2024: Initial coding;
% DJ, 01/03/2025: Display zero equation for x, y, or z as well, but only if
%                   explicitly declared using 'is_zero' field;
% DJ, 11/06/2026: Get the subscripts and variable names from
%                   'varname2display', display the time derivative as d_t
%                   unless 'use_unicode', and call coefficients C_ij;
% MMP, 11/19/2026: Release comment pass (comment-only).
% DJ, 11/23/2026: Give the hint on extracting C_ij for coefficients in any
%                   term or factor, as name.C{i,j} or name.C{i,j,k};

% % % Declare some UNICODE char symbols to display temporal derivatives of
% % % state components for the left-hand side of the equation. The names,
% % % subscripts and superscripts are given by 'varname2display'.
[~,sub_dig,sup_dig] = varname2display('0123456789');
sub_num = cellstr(sub_dig(:));      % Subscripts
sup_num = cellstr(sup_dig(:));      % Superscripts
[~,sub_i] = varname2display('i');
[~,sub_j] = varname2display('j');
sub_t = char(8348);
if ~use_unicode()
    sub_t = '_t';           % the subscript t is missing from many fonts
end
partial = char(8706);
% Name of the variable in the hint on extracting coefficients
if nargin<2 || isempty(name),   name = 'PDE';   end                         % DJ, 11/23/2026


% % % Loop over all types of equations in the PDE, displaying them line by
% % % line.
% Whether coefficients C_ij (single factor) and C_ijk (factor k of a
% product term) are displayed
use_C = false(1,2);                                                         % DJ, 11/23/2026
eq_types = {'free';'x';'y';'z';'BC'};
% If free terms are available, assume there are actual full equations
if numel(PDE.free)>0
    eq_types = {'free'};
end
n_eqs_tot = 0; 
for kk=1:numel(eq_types)
    eq_type = eq_types{kk};    
    if numel(PDE.(eq_type))==0
        % There are no equations of the specified type.
        continue
    end
    % Initialize the left-hand side of the equation.
    if strcmp(eq_type,'x') || strcmp(eq_type,'y') || strcmp(eq_type,'z')
        % For PDE output equations, display an equal between left- and 
        % right-hand side.
        Lobj = eq_type;
        eq_sign = ' = ';
    elseif strcmp(eq_type,'BC')
        % For boundary conditions, display an equal between left- and 
        % right-hand side.
        Lobj = '0';
        eq_sign = ' = ';
        use_ID = 0.5;   % Use the ID only if available.
    else
        % For free terms, there is no equation, so no equal sign.
        Lobj = '';
        eq_sign = '   ';
        use_ID = true;
    end
    % % Loop over all the equations of the considered type.
    for ii=1:numel(PDE.(eq_type))
        eq_struct = PDE.(eq_type){ii};
        if ~isfield(eq_struct,'term') || isempty(eq_struct.term)
            % There are no terms to display in the equation
            if strcmp(eq_type,'free')
                % For free terms, print a zero.
                fprintf('\n    0')
                continue
            elseif ~isfield(eq_struct,'is_zero') || ~eq_struct.is_zero      % DJ, 01/03/2025
                % For equations not explicitly declared as zero equations, 
                % print nothing.
                continue
            end
        end
        
        % % Display the left-hand side of the equation
        if strcmp(eq_type,'x') || strcmp(eq_type,'y') || strcmp(eq_type,'z')
            % % For states, display the temporal derivative
            if ~isfield(eq_struct,'tdiff') || eq_struct.tdiff==0
                Lstate_trm = Lobj;
            elseif eq_struct.tdiff==1
                % First order temporal derivative does not need a superscript.
                Lstate_trm = [partial,sub_t,' ',Lobj];
            elseif eq_struct.tdiff<=9
                % Display a single-digit order derivative.
                Lstate_trm = [partial,sub_t,sup_num{eq_struct.tdiff+1},Lobj];
            else
                % Display a multi-digit order derivative.
                tD_sup = cell2mat(sup_num(str2num(num2str(eq_struct.tdiff)')+1)');
                Lstate_trm = [partial,sub_t,tD_sup,Lobj];
            end
            % % Assign an index to the state or output signal.
            if isfield(eq_struct,'ID')
                Lobj_ID = eq_struct.ID;
                use_ID = true;
            else
                Lobj_ID = ii;
                use_ID = false;
            end
            if Lobj_ID<=9
                Lobj_ID_sub = sub_num{Lobj_ID+1};
            else
                Lobj_ID_sub = cell2mat(sub_num(str2num(num2str(Lobj_ID)')+1)');
            end
            Lstate_trm = [Lstate_trm,Lobj_ID_sub];
            % % Add the variables on which the object depends.
            Lobj_vars = '(t';
            if ~isfield(eq_struct,'vars') || isempty(eq_struct.vars)
                Lobj_vars = [Lobj_vars,')'];
            else
                for ll=1:size(eq_struct.vars,1)
                    Lobj_vars = [Lobj_vars,',',varname2display(eq_struct.vars(ll).varname{1})];
                end
                Lobj_vars = [Lobj_vars,')'];
            end
            Lstate_trm = [Lstate_trm,Lobj_vars];
        else
            Lstate_trm = Lobj;
        end
        % % Display the left-hand side of the equation.
        fprintf(['\n ',Lstate_trm,eq_sign])
    
        % % Print the first term
        if ~isfield(eq_struct,'term') || isempty(eq_struct.term)            % DJ, 01/03/2025
            % If there is no term, the RHS is just zero.
            fprintf('0;');
            continue
        end
        % Extract the first factor of the term
        if isscalar(eq_struct.term{1})
            C_idx = [ii+n_eqs_tot,1];
        else
            C_idx = [ii+n_eqs_tot,1,1];
        end
        % The product sign: '*' for a scalar equation, '.*' otherwise
        if isfield(eq_struct,'size') && eq_struct.size==1
            prod_str = '*';
        else
            prod_str = '.*';
        end
        % The factors of the first term joined by the product sign, the
        % sign of the term, and whether any factor uses C_ij
        [term_str,term_sign,use_C1] = display_term(PDE,eq_struct.term{1}(1),C_idx,use_ID);    % DJ, 11/23/2026
        for k=2:numel(eq_struct.term{1})
            [term_str_k,term_sign_k,use_Ck] = display_term(PDE,eq_struct.term{1}(k),[ii+n_eqs_tot,1,k],use_ID);
            term_sign = term_sign*term_sign_k;
            term_str = [term_str,prod_str,term_str_k];
            use_C1 = use_C1 || use_Ck;                                      % DJ, 11/23/2026
        end
        % C_idx has 2 indices (C_ij) or 3 (C_ijk)
        use_C(numel(C_idx)-1) = use_C(numel(C_idx)-1) || use_C1;            % DJ, 11/23/2026
        if term_sign==-1
            fprintf(['- ',term_str]);
        else
            fprintf(term_str);
        end
        % % Print the remaining terms.
        for jj=2:numel(eq_struct.term)
            % Extract the first factor of the term
            if isscalar(eq_struct.term{jj})
                C_idx = [ii+n_eqs_tot,jj];
            else
                C_idx = [ii+n_eqs_tot,jj,1];
            end
            % The factors of term jj, its sign, and whether it uses C_ij
            [term_str,term_sign,use_Cjj] = display_term(PDE,eq_struct.term{jj}(1),C_idx,use_ID);
            for k=2:numel(eq_struct.term{jj})
                [term_str_k,term_sign_k,use_Ck] = display_term(PDE,eq_struct.term{jj}(k),[ii+n_eqs_tot,jj,k],use_ID);
                term_sign = term_sign*term_sign_k;
                term_str = [term_str,prod_str,term_str_k];
                use_Cjj = use_Cjj || use_Ck;                                % DJ, 11/23/2026
            end
            if term_sign==-1
                fprintf([' - ',term_str]);
            else
                fprintf([' + ',term_str]);
            end
            use_C(numel(C_idx)-1) = use_C(numel(C_idx)-1) || use_Cjj;       % DJ, 11/23/2026
        end
        fprintf(';');
    end
    n_eqs_tot = n_eqs_tot + numel(PDE.(eq_type));
end

% % Inform users on how to extract coefficients, if applicable, with the
% % calls for the numbers of indices displayed.
if any(use_C)                                                               % DJ, 11/23/2026
    [~,sub_k] = varname2display('k');
    calls = {['"',name,'.C{i,j}"'],['"',name,'.C{i,j,k}"']};
    coeffs = {['C',sub_i,sub_j],['C',sub_i,sub_j,sub_k]};
    fprintf('\n\n Call %s to see the value of coefficients %s as displayed.',strjoin(calls(use_C),' or '),strjoin(coeffs(use_C),', '));
end

% % Add some empty lines.
fprintf('\n\n')

end

