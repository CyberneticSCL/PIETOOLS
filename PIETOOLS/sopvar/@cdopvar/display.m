function display(Pop,name,width_max,expand)
% DISPLAY(POP,NAME,WIDTH_MAX,EXPAND) displays the block matrix POP of PI
% operators with decision variables in the command window as
%   NAME = [C_11 ... C_1d] : [X_1] --> [Y_1]
%          [ :   .    :  ]   [ : ]     [ : ]
%          [C_p1 ... C_pd]   [X_d]     [Y_p]
% with each of the blocks C_ij named, and zero blocks displayed as 0, and
% with X_j the space that column j acts on, and Y_i the space that row i
% maps to, e.g. R, L_2[0,1], or L_2^2[[-1,1]x[0,2]], with the 2 of L_2 as
% subscript and the powers as superscripts, followed by the operator
% written out explicitly, or the definitions of the blocks as the actual
% operators acting on x_j, if they are small enough, see 'get_display'. If
% the operator maps from a single space to a single space, it is displayed
% as
%   NAME : X --> Y
% instead, followed by the single block, see 'sdopvar/display'. An empty
% block matrix is displayed as
%   NAME = [] : [X_1] --> [Y_1]
% with R^0 as the space of the absent block columns or rows. Called
% automatically when displaying a 'cdopvar' object without semicolon.
%
% INPUTS
% - Pop:    'cdopvar' object;
% - name:   (optional) char, the name to display for the operator. Defaults
%           to the name of the variable passed to display, or to 'Pop' if
%           there is none. The name 'ans' is replaced by 'Pop', but is
%           displayed as ans = Pop = [C_11 ...] : ... or as
%           ans = Pop : X --> Y;
% - width_max:  (optional) maximal number of characters of a line with
%           multiple terms. Defaults to the width of the command window, or
%           to 100 if this cannot be determined;
% - expand: (optional) logical, false to not define the blocks as the
%           actual operators. Defaults to true;
%
% OUTPUTS
% Display of the block matrix in the command window, with the lines
% constructed by 'get_display'.
%
% NOTES
% - An identity operator is displayed as (NAME x)_i(s) = x_i(s) for each
%   block row i, with the spatial variables of the block, and a block that
%   is an identity operator is defined as (C_ij x_j)(s) = x_j(s).
% - The blocks are named C_ij, dropping i if there is only one block row,
%   and j if there is only one block column.
% - The blocks can be extracted by calling Pop.C{i,j}, and displayed in
%   turn. The decision variables of the blocks are those in Pop.Zd.
% - The spaces follow from Pop.dim_in, Pop.dim_out, Pop.space_in,
%   Pop.space_out, Pop.vars and Pop.dom, see 'get_display'.

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
% DJ, 11/06/2026: Initial coding;

% Determine the name of the operator
if nargin<2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'Pop';
    end
end
call_name = name;                       % the name that the user can call
is_ans = strcmp(name,'ans');
if is_ans
    % Display an unnamed object under a default name
    name = 'Pop';
end
if nargin<3
    width_max = [];
end
if nargin<4
    expand = true;
end
[p,d] = size(Pop.C);

% Construct the lines of the display, and the spaces. A single block that
% is displayed by its own display function is not defined here
is_single = p==1 && d==1;
is_deleg = is_single && (isa(Pop.C{1,1},'sdopvar') || isa(Pop.C{1,1},'sopvar'));
disp_struct = get_display(Pop,name,expand && ~is_deleg);

% A single block, between single spaces, is displayed as NAME : X --> Y,
% followed by the block
if is_single
    ans_str = '';
    if is_ans
        ans_str = 'ans = ';
    end
    fprintf('\n  %s%s : %s --> %s\n',ans_str,name,disp_struct.spc_in{1},disp_struct.spc_out{1});
end
if is_deleg
    display(Pop.C{1,1},name,width_max,[call_name,'.C{1,1}']);
    return
end

% Display the matrix of blocks next to the spaces
if disp_struct.is_id
    % An identity operator is displayed as (Pop x)_i(s) = x_i(s) for each row i
    if is_ans && ~is_single
        fprintf('\n  ans = %s, where\n',name);
    end
    fprintf('\n');
    for i=1:p
        args = '';
        if ~isempty(Pop.C{i,i}.vars.out)
            out_names = cellfun(@varname2display,reshape(Pop.C{i,i}.vars.out,1,[]),'UniformOutput',false);
            args = ['(',strjoin(out_names,','),')'];
        end
        lhs = [name,' ',disp_struct.xsym];
        if p>1
            lhs = ['(',lhs,')',sprintf('%d',i)-'0'+8320];
        end
        print_terms([lhs,args,' ='],{[disp_struct.xnames{i},args]},1,width_max);
    end
    fprintf('\n');
    return
end
if ~is_single
    % Merge ans = into the first line, and indent the others to match
    lines = disp_struct.lines;
    if is_ans
        lines{1} = ['  ans = ',lines{1}(3:end)];
        lines(2:end) = cellfun(@(a) ['      ',a],lines(2:end),'UniformOutput',false);
    end
    fprintf('\n');
    fprintf('%s\n',lines{:});
end
fprintf('\n');
if p==0 || d==0
    % No blocks to define or extract
    return
end

% Display the operator explicitly, or the definitions of the blocks
for k=1:numel(disp_struct.rows)
    print_terms(disp_struct.rows(k).lhs,disp_struct.rows(k).terms,disp_struct.rows(k).signs,width_max);
end
if ~isempty(disp_struct.rows)
    fprintf('\n');
end
for k=1:numel(disp_struct.defs)
    print_terms(disp_struct.defs(k).lhs,disp_struct.defs(k).terms,disp_struct.defs(k).signs,width_max);
end
if ~isempty(disp_struct.defs)
    fprintf('\n');
end

% Explain how to extract the blocks
blk_name = 'C';
if p>1
    blk_name = [blk_name,char(7522)];       % subscript i
end
if d>1
    blk_name = [blk_name,char(11388)];      % subscript j
end
if p>1
    iname = 'i';
else
    iname = '1';
end
if d>1
    jname = 'j';
else
    jname = '1';
end
fprintf('  Call "%s.C{%s,%s}" to extract block %s.\n\n',call_name,iname,jname,blk_name);

end
