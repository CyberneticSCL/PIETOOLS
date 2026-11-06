function display(PIE,name,width_max,expand)
% DISPLAY(PIE,NAME,WIDTH_MAX,EXPAND) displays the equations
% defining the PIE PIE in the command window, as
%   dt(T x^ + Tw w^ + Tu u^)(t) = A x^(t) + Bw w^(t) + Bu u^(t)
%                          z(t) = Cz x^(t) + Dzw w^(t) + Dzu u^(t)
%                          y(t) = Cy x^(t) + Dyw w^(t) + Dyu u^(t)
%                             0 = K x^(t) + Kw w^(t) + Ku u^(t),
% with dt the partial derivative in time, x^, w^, u^ the fundamental state
% and inputs, and the operators defined as the actual operators acting on
% the state and inputs if they are small enough. For a nonlinear PIE
% ('npie_struct'), A x^, Cz x^ and Cy x^ are replaced by the distributed
% polynomials f(x^), gz(x^) and gy(x^).
%
% INPUTS
% - PIE:    'pie_struct' or 'npie_struct' object;
% - name:   (optional) char, the name to display for the PIE. Defaults to
%           the name of the variable passed to display, or to 'PIE' if
%           there is none. The name 'ans' is replaced by 'PIE';
% - width_max:  (optional) maximal number of characters of a line with
%           multiple terms. Defaults to the width of the command window, or
%           to 100 if this cannot be determined;
% - expand: (optional) logical, false to not define the operators as the
%           actual operators. Defaults to true;
%
% OUTPUTS
% Display of the equations in the command window.
%
% NOTES
% - Terms that are zero or empty are omitted. The regulated output z and
%   the observed output y are only displayed if the PIE has any, and the
%   algebraic constraint only if it is nontrivial.
% - Below the equations, after 'where' and indented, the operators are
%   defined as the actual operators acting on the state and inputs, and
%   after an empty line, the space of each nonempty state, input and output
%   is stated, e.g. x^(t) in L_2^3[[0,1]x[-1,1]].
% - The operators are named as the fields of PIE, and defined as the actual
%   operators acting on the state and inputs, using 'get_display' of the
%   classes of the operators, if they are scalar-valued and small enough.
%   The other operators can be extracted by calling PIE.T, PIE.Bw, etc., as
%   indicated below the equations. The terms are placed on the same line,
%   unless the line would have more characters than fit in the command
%   window, in which case the whole next term is moved to the next line,
%   see 'print_terms'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - display
%
% Copyright (C) 2026 PIETOOLS Team
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
% Initial coding DJ - 08/30/2022
% DJ, 06/25/2025: Add display of 'misc' field;
% DB, 10/27/2026: Display optional PI operators F/K in case the null
% space of the derivative operator is non-trivial.
% DJ, 10/30/2026: Update to new pie_struct structure (display the function
%                   spaces, and constraint operators if specified);
% DJ, 11/06/2026: Define the operators as actual operators, if small, and
%                   display the spaces of the state, inputs and outputs;
% DJ, 11/06/2026: Get the names of PIE states and the subscripts from
%                   'varname2display';

if nargin<2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'PIE';
    end
end
call_name = name;
if nargin<3
    width_max = [];
end
if nargin<4
    expand = true;
end
if ~(isa(PIE,'pie_struct') || isa(PIE,'npie_struct'))
    error("The PIE must be specified as 'pie_struct' or 'npie_struct'.")
end
if strcmp(name,'ans')
    % Display an unnamed PIE under a default name
    name = 'PIE';
    fprintf('\n  ans = %s, where\n',name);
end

% An empty PIE, for which initialize fails, has no equations to display
if ~PIE.is_initialized
    try
        PIE = initialize(PIE);
    catch
        fprintf('\n  %s is an empty or incomplete PIE.\n\n',name);
        return
    end
end

% Symbols: the partial derivative in time, and the names of the state, inputs and
% outputs of the PIE, as given by 'varname2display' for a name ending with _pie.
% The inputs w, u are named as for the PIE, unless the map PIE.pie2pde.w,
% PIE.pie2pde.u to the PDE input is specified and an identity operator, so that
% it is the PDE input. The outputs z, y are named as for the PDE, unless the map
% to the PDE output is specified, and not an identity operator
dt_sym = [char(8706),char(8348)];
sig_names = struct('w',varname2display('w_pie'),'u',varname2display('u_pie'),'y','y','z','z');
for sig = {'w','u','y','z'}
    is_specified = false;       is_pde_id = false;
    try
        is_specified = prod(size(PIE.pie2pde.(sig{1})))>0;       % an empty cell, or 0 x 0 tensopmat, if not
        is_pde_id = is_specified && is_id_op(PIE.pie2pde.(sig{1}));
    catch
    end
    if is_pde_id && any(strcmp(sig{1},{'w','u'}))
        sig_names.(sig{1}) = sig{1};
    elseif is_specified && ~is_pde_id && any(strcmp(sig{1},{'y','z'}))
        sig_names.(sig{1}) = varname2display([sig{1},'_pie']);
    end
end
x_names = {varname2display('x_pie'),sig_names.w,sig_names.u};
is_nl = isa(PIE,'npie_struct');

% Operators in the left-hand side, dynamics, regulated output, observed output
% and constraint (rows), acting on the state, exogenous inputs and actuator
% inputs (columns). The polynomials replace A, Cz, Cy for a nonlinear PIE.
eq_ops = {'T','Tw','Tu'; 'A','Bw','Bu'; 'Cz','Dzw','Dzu'; 'Cy','Dyw','Dyu'; 'K','Kw','Ku'};
if is_nl
    eq_ops(2:4,1) = {'f';'gz';'gy'};
end

% Terms of each equation, for the operators that are nonzero
eq_terms = cell(5,1);
is_used = false(5,3);
for r=1:5
    eq_terms{r} = cell(1,0);
    for c=1:3
        if is_zero(PIE.(eq_ops{r,c}))
            continue
        end
        is_used(r,c) = true;
        if r==1
            eq_terms{r}{end+1} = [eq_ops{r,c},' ',x_names{c}];
        elseif is_nl && c==1 && r<=4
            eq_terms{r}{end+1} = [eq_ops{r,c},'(',x_names{1},'(t))'];
        else
            eq_terms{r}{end+1} = [eq_ops{r,c},' ',x_names{c},'(t)'];
        end
    end
end
if isempty(eq_terms{1})
    eq_terms{1} = {'0'};
end

% Display the equations with the equality signs aligned. The outputs are
% displayed if there are any, and the constraint if it is nontrivial. The width
% of a left-hand side excludes the combining marks, which take no space.
eq_lhs = {[dt_sym,' (',strjoin(eq_terms{1},' + '),')(t) ='], [sig_names.z,'(t) ='], [sig_names.y,'(t) ='], '0 ='};
eq_shown = [true, ~isempty(PIE.space.z.dims), ~isempty(PIE.space.y.dims), ~isempty(eq_terms{5})];
lhs_widths = cellfun(@(a) numel(a)-nnz(a>=768 & a<=879),eq_lhs);     % combining marks U+0300-U+036F
lhs_widths(1) = lhs_widths(1)-1;        % the space after dt_sym makes up for the narrow subscript t
pie_type = 'linear';
if is_nl
    pie_type = 'nonlinear';
end
fprintf('\n  %s is a %s PIE of the form\n\n',name,pie_type);
for q=find(eq_shown)
    print_terms([repmat(' ',1,max(lhs_widths(eq_shown))-lhs_widths(q)),eq_lhs{q}],eq_terms{q+1},[],width_max);
end
fprintf('\n');

% The spaces of the state, inputs and outputs that are nonempty: the product
% of the spaces L_2^m[S] of the components, with S the product of the
% intervals of the variables that they depend on, and R^m if there are none
signals = {'x',[x_names{1},'(t)']; 'w',[x_names{2},'(t)']; 'u',[x_names{3},'(t)']; 'z',[sig_names.z,'(t)']; 'y',[sig_names.y,'(t)']};
space_lines = cell(1,0);
for k=1:size(signals,1)
    sig_space = PIE.space.(signals{k,1});
    comp_strs = cell(1,numel(sig_space.dims));
    for l=1:numel(sig_space.dims)
        pow_str = '';
        if sig_space.dims(l)~=1
            [~,~,pow_str] = varname2display(sprintf('%d',sig_space.dims(l)));       % superscript
        end
        var_idcs = find(sig_space.depmat(l,:));
        if isempty(var_idcs)
            comp_strs{l} = ['R',pow_str];
        else
            intvl_strs = arrayfun(@(q) sprintf('[%g,%g]',PIE.space.dom(q,1),PIE.space.dom(q,2)),var_idcs,'UniformOutput',false);
            if ~isscalar(var_idcs)
                intvl_strs = {['[',strjoin(intvl_strs,'x'),']']};
            end
            comp_strs{l} = [varname2display('L2'),pow_str,intvl_strs{1}];
        end
    end
    if ~isempty(comp_strs)
        space_lines{end+1} = [signals{k,2},' in ',strjoin(comp_strs,[' ',char(215),' '])];      %#ok<AGROW>
    end
end

% Define the operators as the actual operators acting on the state and inputs,
% if they are scalar-valued and small enough, and the polynomials as the sums
% of operators C_j acting on monomials
def_lhs = cell(1,0);        def_terms = cell(1,0);      def_signs = cell(1,0);
hidden = cell(1,0);                         % operators that are not defined
poly_names = cell(1,0);                     % polynomials with named operators C_j
poly_nterms = zeros(1,0);
[used_c,used_r] = find(is_used.');
for q=1:numel(used_r)
    op_name = eq_ops{used_r(q),used_c(q)};
    x_name = x_names{used_c(q)};
    pie_op = PIE.(op_name);
    is_defined = false;
    is_named = false;
    out_args = cell(1,0);
    try
        if isa(pie_op,'polyopvar')
            % A polynomial with a single row, in the state x^ if there is only one
            is_defined = size(pie_op.C.ops,1)==1;
            if is_defined
                if ~pie_op.chkval
                    pie_op = initialize(pie_op);
                end
                if isscalar(PIE.space.x.dims) && isscalar(pie_op.varname)
                    pie_op.varname = {'x_pie'};
                end
                blk_disp = get_display(pie_op,expand);
                x_name = '';
                if ~isempty(blk_disp.states)
                    x_name = ['(',strjoin(blk_disp.states,','),')'];
                end
                out_args = blk_disp.out_args;
                is_named = ~blk_disp.expanded;
            end
        elseif expand && isa(pie_op,'copvar') && isequal(size(pie_op.C),[1,1]) && ...
                isequal(pie_op.dim_in,1) && isequal(pie_op.dim_out,1) && isa(pie_op.C{1,1},'sopvar')
            blk_disp = get_display(pie_op.C{1,1},x_name);
            is_defined = blk_disp.ok;
            out_args = cellfun(@varname2display,reshape(pie_op.C{1,1}.vars.out,1,[]),'UniformOutput',false);
        elseif expand && isa(pie_op,'tensopmat') && isequal(size(pie_op.ops),[1,1]) && ...
                isequal(pie_op.dims.in,1) && isequal(pie_op.dims.out,1) && ...
                isa(pie_op.ops{1,1},'tensopvar') && size(pie_op.ops{1,1}.ops,2)==1
            blk_disp = get_display(pie_op.ops{1,1},{x_name});
            is_defined = blk_disp.expanded;
            out_args = blk_disp.out_args;
        end
        if is_defined
            % Left-hand side (NAME x)(s) = ..., or NAME(x)(s) = ... for a polynomial
            lhs_str = [op_name,' ',x_name];
            if isa(pie_op,'polyopvar')
                lhs_str = [op_name,x_name];
            elseif ~isempty(out_args)
                lhs_str = ['(',lhs_str,')'];
            end
            if ~isempty(out_args)
                lhs_str = [lhs_str,'(',strjoin(out_args,','),')'];
            end
            def_lhs{end+1} = [lhs_str,' ='];                %#ok<AGROW>
            def_terms{end+1} = blk_disp.terms;              %#ok<AGROW>
            def_signs{end+1} = blk_disp.signs;              %#ok<AGROW>
            if is_named
                poly_names{end+1} = op_name;                %#ok<AGROW>
                poly_nterms(end+1) = numel(blk_disp.terms); %#ok<AGROW>
            end
        end
    catch
        is_defined = false;
    end
    if ~is_defined
        hidden{end+1} = op_name;                            %#ok<AGROW>
    end
end

% Define the operators, and state the spaces after an empty line, both indented
% below 'where'
if ~isempty(space_lines) || ~isempty(def_lhs)
    fprintf('  where\n');
    for q=1:numel(def_lhs)
        print_terms(['  ',def_lhs{q}],def_terms{q},def_signs{q},width_max);
    end
    if ~isempty(def_lhs) && ~isempty(space_lines)
        fprintf('\n');
    end
    if ~isempty(space_lines)
        fprintf('    %s\n',space_lines{:});
    end
    fprintf('\n');
end

% Explain how to extract the operators C_j of the polynomials that are not
% displayed, and the operators that are not defined
for q=1:numel(poly_names)
    jname = '1';
    op_name = 'C';
    if poly_nterms(q)>1
        jname = 'j';
        [~,sub_j] = varname2display('j');
        op_name = ['C',sub_j];
    end
    fprintf('  Call "%s.%s.C{1,%s}" to extract operator %s acting on monomial %s of %s.\n', ...
        call_name,poly_names{q},jname,op_name,jname,poly_names{q});
end
if ~isempty(hidden)
    calls = strcat('"',call_name,'.',hidden,'"');
    if numel(hidden)>3
        fprintf('  Call %s, etc., to extract the different operators.\n',strjoin(calls(1:3),', '));
    else
        fprintf('  Call %s to extract %s.\n',strjoin(calls,', '),strjoin(hidden,', '));
    end
end
if ~isempty(poly_names) || ~isempty(hidden)
    fprintf('\n');
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function is_zero_op = is_zero(obj)
% IS_ZERO_OP = IS_ZERO(OBJ) returns true if the operator OBJ of a PIE, a
% 'copvar', 'tensopmat' or 'polyopvar' object, or a block of such an
% operator, is empty or zero.
%
% INPUTS
% - obj:    'copvar', 'tensopmat', 'polyopvar', 'tensopvar', 'sopvar' or
%           'intop' object, or an array of doubles;
%
% OUTPUTS
% - is_zero_op: logical, true if all blocks of obj are empty or zero. A
%           'tensopvar' is zero if every term has a zero factor;
if isa(obj,'copvar')
    is_zero_op = all(cellfun(@is_zero,obj.C(:)));
elseif isa(obj,'tensopmat')
    is_zero_op = all(cellfun(@is_zero,obj.ops(:)));
elseif isa(obj,'polyopvar')
    is_zero_op = all(cellfun(@is_zero,obj.C.ops(:)));
elseif isa(obj,'tensopvar')
    is_zero_op = all(any(cellfun(@is_zero,obj.ops),2));
elseif isa(obj,'sopvar') || isa(obj,'intop')
    is_zero_op = all(cellfun(@nnz,obj.params(:))==0);
else
    is_zero_op = isempty(obj) || ~any(obj(:));
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %