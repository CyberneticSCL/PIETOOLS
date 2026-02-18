function [] = display_PDE_crude(PDE)
% DISPLAY_PDE_CRUDE(PDE) displays the input PDE in the Command Window,
% giving a symbolic overview of the specified equations/terms using
% UNICODE. This function is intended to display "unfinished"
% (un-initialized) pde_struct objects, as the standard "display_PDE"
% function may display initialized structures more cleanly.
%
% INPUT
% - PDE:    'pde_struct' object, representing a PDE to be displayed.
%           Equations will be printed in the Command Window.
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
% DJ, 06/23/2024: Initial coding;
% DJ, 01/03/2025: Display zero equation for x, y, or z as well, but only if
%                   explicitly declared using 'is_zero' field;

% % % Declare some UNICODE char symbols to display temporal derivatives of
% % % state components for the left-hand side of the equation.
% % https://rupertshepherd.info/resource_pages/superscript-letters-in-unicode
sub_i = '\x1D62';
%sub_j = '\x2C7C';
sub_k = '\x2096';
sub_t = '\x209C';
partial = '\x2202';
sub_num = mat2cell([repmat('\x208',[10,1]),num2str((0:9)')],ones(10,1),6);
sup_num = cell(10,1);    % Superscripts
sup_num{1} = '\x2070';   % Superscript 0
sup_num{2} = '\xB9';
sup_num{3} = '\xB2';
sup_num{4} = '\xB3';
sup_num{5} = '\x2074';
sup_num{6} = '\x2075';
sup_num{7} = '\x2076';
sup_num{8} = '\x2077';
sup_num{9} = '\x2078';
sup_num{10} = '\x2079';
sup_min = '\x207B';


% % % Loop over all types of equations in the PDE, displaying them line by
% % % line.
use_C = false;
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
                    Lobj_vars = [Lobj_vars,',',eq_struct.vars(ll).varname{1}];
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
        if isfield(eq_struct,'size') && eq_struct.size==1
            prod_str = '*';
        else
            prod_str = '.*';
        end
        [term_str,term_sign] = display_term(PDE,eq_struct.term{1}(1),C_idx,use_ID);
        for k=2:numel(eq_struct.term{1})
            [term_str_k,term_sign_k] = display_term(PDE,eq_struct.term{1}(k),[ii+n_eqs_tot,1,k],use_ID);
            term_sign = term_sign*term_sign_k;
            term_str = [term_str,prod_str,term_str_k];
        end
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
            [term_str,term_sign,use_Cjj] = display_term(PDE,eq_struct.term{jj}(1),C_idx,use_ID);
            for k=2:numel(eq_struct.term{jj})
                [term_str_k,term_sign_k] = display_term(PDE,eq_struct.term{jj}(k),[ii+n_eqs_tot,jj,k],use_ID);
                term_sign = term_sign*term_sign_k;
                term_str = [term_str,prod_str,term_str_k];
            end
            if term_sign==-1
                fprintf([' - ',term_str]);
            else
                fprintf([' + ',term_str]);
            end
            use_C = use_C || use_Cjj;
        end
        fprintf(';');
    end
    n_eqs_tot = n_eqs_tot + numel(PDE.(eq_type));
end

% % Inform users on how to extract coefficients, if applicable.
if use_C
    fprintf(['\n\n Call "PDE.C{i,k}" to see the value of coefficients C',sub_i,sub_k,' as in the displayed equations.'])
end

% % Add some empty lines.
fprintf('\n\n')

end

