function disp_struct = get_display(Pop,name,expand)
% DISP_STRUCT = GET_DISPLAY(POP,NAME,EXPAND) constructs the lines of the
% display of the block matrix POP of PI operators, for use in the 'display'
% method of 'copvar',
%   NAME = [C_11 ... C_1d] : [X_1] --> [Y_1]
%          [ :   .    :  ]   [ : ]     [ : ]
%          [C_p1 ... C_pd]   [X_d]     [Y_p]
% with X_j the space that column j acts on, and Y_i the space that row i
% maps to, and the operator written out, or the definitions of the blocks,
% if they are small enough.
%
% INPUTS
% - Pop:    'copvar' object, with p x d blocks;
% - name:   (optional) char, the name of the block matrix. Defaults to
%           'Pop';
% - expand: (optional) logical, true if the blocks may be defined as the
%           actual operators if they are small enough. Defaults to true;
%
% OUTPUTS
% - disp_struct:      struct with the fields
%       lines:  cellstr, the lines of the display of the matrix of blocks
%               next to the spaces;
%       spc_in: cellstr, the spaces that the block columns act on, e.g.
%               'L_2[0,1]';
%       spc_out:    cellstr, the spaces that the block rows map to;
%       defs:   struct array with the fields lhs, terms and signs, for each
%               of the blocks that is displayed as the actual operator, in
%               the order of the rows, with lhs the left-hand side, e.g.
%               '(C_11 x_1)(s) =', and terms and signs as for the output D
%               of 'sopvar/get_display'. An identity block is defined with
%               its spatial variables, (C_11 x_1)(s) = x_1(s);
%       rows:   struct array with the fields lhs, terms and signs, for each
%               element of the output of Pop, with lhs the left-hand side,
%               e.g. '(Pop x)_1(s) =', and terms and signs for the sum of
%               the blocks in its row acting on the elements of the input,
%               see 'get_display'. Empty if the operator is not displayed
%               explicitly, see the notes;
%       xnames: 1 x d cellstr, the names x_j of the functions that the
%               blocks act on;
%       xsym:   char, the name x of the function that Pop acts on, from
%               which x_j follows;
%       is_id:  logical, true if Pop is an identity operator, in which case
%               'display' shows (Pop x)_i(s) = x_i(s) for each block
%               row i;
%
% NOTES
% - The blocks are named C_ij, dropping i if there is only one block row,
%   and j if there is only one block column. A zero block is displayed as
%   0.
% - The spaces follow from Pop.dim_in, Pop.dim_out, Pop.space_in,
%   Pop.space_out, Pop.vars and Pop.dom. A space is L_2^m[S] for m the
%   dimension (dropped for m=1) and S the domain of the functions, the
%   product of the intervals of the variables that they depend on. With a
%   single interval, S is that interval, e.g. L_2[0,1]. If the functions
%   depend on no variables, the space is R^m, or R for m=1. The names of
%   the variables are not displayed. The subscript 2 and the powers are
%   displayed as subscript and superscript characters.
% - An empty block matrix (p=0 or d=0) is displayed as NAME = [] between
%   the spaces, with the space R^0 for the absent block columns or rows,
%   and nothing is defined.
% - The operator is displayed explicitly, with a row for each element of
%   the output, if the input and the output have at most 5 elements,
%   counting those of vector-valued functions, every block is displayed
%   explicitly by 'get_display' (scalar-valued, or with explicit rows), and
%   each row has at most 6 terms of at most 100 characters. The elements of
%   a function x_j of a vector-valued space are named x_j,k. The blocks are
%   then not defined.
% - The definitions of the blocks are only given if all of the blocks can
%   be displayed as the actual operators, acting on x_j, with 'get_display'
%   of 'sopvar', and there are at most 6 of them. Blocks that are
%   matrix-valued are not displayed as the actual operators, so that none
%   of the blocks is defined if one of them is.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_display
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
% DJ, 11/06/2026: Get the subscripts from 'varname2display';

% Default settings
max_blocks = 6;         % maximal number of blocks defined as the actual operators
max_elems = 5;          % maximal number of elements of the input and of the output
                        % for which the operator is displayed explicitly
max_terms = 6;          % maximal number of terms of each element of the output
max_chars = 100;        % maximal number of characters of each such term

if nargin<2
    name = 'Pop';
end
if nargin<3
    expand = true;
end
[p,d] = size(Pop.C);

[~,sub_dig,sup_dig] = varname2display('0123456789');        % subscript and superscript digits
sub_str = @(n) sub_dig(sprintf('%d',n)-'0'+1);
sup_str = @(n) sup_dig(sprintf('%d',n)-'0'+1);
L2_sym = ['L',sub_str(2)];

blk = cell(p,d);
for i=1:p
    for j=1:d
        if isempty(Pop.C{i,j})
            blk{i,j} = '0';
            continue
        end
        blk_name = 'C';
        if p>1
            blk_name = [blk_name,sub_str(i)];
        end
        if p>1 && d>1 && max(i,j)>9
            blk_name = [blk_name,','];
        end
        if d>1
            blk_name = [blk_name,sub_str(j)];
        end
        blk{i,j} = blk_name;
    end
end

% Determine the spaces that the block columns act on (s=1) and that the block
% rows map to (s=2), from the dimensions and the variables in each space
dims = {Pop.dim_in,Pop.dim_out};
masks = {Pop.space_in,Pop.space_out};
spc = cell(1,2);
for s=1:2
    spc{s} = cell(numel(dims{s}),1);
    for l=1:numel(dims{s})
        m = dims{s}(l);
        var_idcs = find(masks{s}(l,:));
        if isempty(var_idcs)
            % Functions of no variables: constants
            spc_str = 'R';
            if m~=1
                spc_str = ['R',sup_str(m)];
            end
        else
            % Product of the intervals of the variables
            dom_strs = cell(1,numel(var_idcs));
            for q=1:numel(var_idcs)
                dom_strs{q} = sprintf('[%g,%g]',Pop.dom(var_idcs(q),1),Pop.dom(var_idcs(q),2));
            end
            if isscalar(var_idcs)
                dom_str = dom_strs{1};
            else
                dom_str = ['[',strjoin(dom_strs,'x'),']'];
            end
            spc_str = [L2_sym,dom_str];
            if m~=1
                spc_str = [L2_sym,sup_str(m),dom_str];
            end
        end
        spc{s}{l} = spc_str;
    end
end

% An empty block matrix is displayed as [], with the space R^0 for the absent
% block columns or rows
is_empty = p==0 || d==0;
for s=1:2
    if is_empty && isempty(spc{s})
        spc{s} = {['R',sup_str(0)]};
    end
end

% Lines of the display of the block matrix next to the spaces that it maps from
% and to,
%   NAME = [C_11 ... C_1d] : [X_1] --> [Y_1]
%          [ :   .    :  ]   [ : ]     [ : ]
%          [C_p1 ... C_pd]   [X_d]     [Y_p]
% as max(p,d) lines without trailing white space. The names of the blocks
% are padded to a common width, and so are the spaces, each displayed between
% square brackets as an entry of a column vector. The arrow is only on the
% first line
nrows = max(numel(spc{1}),numel(spc{2}));

blk_width = max(cellfun(@numel,blk(:)));
row_str = cell(p,1);
for i=1:p
    row_str{i} = ['[ ',strjoin(cellfun(@(a) [a,repmat(' ',1,blk_width-numel(a))],blk(i,:),'UniformOutput',false),'  '),' ]'];
end
if is_empty
    row_str = {'[]'};
end
spc_width = [max(cellfun(@numel,spc{1})),max(cellfun(@numel,spc{2}))];
spc_in_str = cellfun(@(a) ['[ ',a,repmat(' ',1,spc_width(1)-numel(a)),' ]'],spc{1},'UniformOutput',false);
spc_out_str = cellfun(@(a) ['[ ',a,repmat(' ',1,spc_width(2)-numel(a)),' ]'],spc{2},'UniformOutput',false);

% The row of blocks, the input space of each block column, an arrow, and the
% output space of each block row, as column vectors
lhs = [name,' = '];
lines = cell(nrows,1);
for k=1:nrows
    if k==1
        line_str = ['  ',lhs];
    else
        line_str = repmat(' ',1,2+numel(lhs));
    end
    if k<=numel(row_str)
        line_str = [line_str,row_str{k}];
    else
        line_str = [line_str,repmat(' ',1,numel(row_str{1}))];
    end
    if k<=numel(spc{1})
        spc_in_k = spc_in_str{k};
    else
        spc_in_k = repmat(' ',1,spc_width(1)+4);
    end
    spc_out_k = '';
    if k<=numel(spc{2})
        spc_out_k = spc_out_str{k};
    end
    if k==1
        line_str = [line_str,' : ',spc_in_k,' --> ',spc_out_k];
    else
        line_str = [line_str,'   ',spc_in_k,'     ',spc_out_k];
    end
    lines{k} = regexprep(line_str,'\s+$','');
end

disp_struct = struct('lines',{lines}, ...
    'defs',struct('lhs',{},'terms',{},'signs',{}),'xnames',{cell(1,d)},'xsym','','is_id',false, ...
    'rows',struct('lhs',{},'terms',{},'signs',{}),'spc_in',{spc{1}},'spc_out',{spc{2}});

% Name of the functions on which the blocks act, x, unless this is also the
% name of a variable (ignoring a number at the end of the name), in which case
% the first of u, w, v that is not
var_bases = regexprep(reshape(Pop.vars,1,[]),'_?\d+$','');
for x_cand = 'xuwv'
    xsym = x_cand;
    if ~ismember(x_cand,var_bases)
        break
    end
end
for j=1:d
    disp_struct.xnames{j} = xsym;
    if d>1
        disp_struct.xnames{j} = [xsym,sub_str(j)];
    end
end
disp_struct.xsym = xsym;
disp_struct.is_id = is_id_op(Pop);
if is_empty
    % Nothing to define
    return
end

% Display the operator explicitly, with each element of the output as the sum
% of the blocks acting on the elements of the input, if there are few
% elements, and each element has few terms of few characters. The elements of
% a vector-valued function x_j are named x_j,k.
if expand && sum(Pop.dim_in)<=max_elems && sum(Pop.dim_out)<=max_elems
    rows = struct('lhs',{},'terms',{},'signs',{});
    is_explicit = true;
    for i=1:p
        % Names of the output variables of the row, from its first nonempty block
        out_names = cell(1,0);
        j = find(~cellfun(@isempty,Pop.C(i,:)),1);
        if ~isempty(j)
            out_names = cellfun(@varname2display,reshape(Pop.C{i,j}.vars.out,1,[]),'UniformOutput',false);
        end
        for a=1:Pop.dim_out(i)
            terms = cell(1,0);      signs = zeros(1,0);
            for j=1:d
                if isempty(Pop.C{i,j})
                    continue
                end
                blk_disp = get_display(Pop.C{i,j},disp_struct.xnames{j});
                if blk_disp.is_id
                    % Identity operator: the element of the function itself
                    x_elem = disp_struct.xnames{j};
                    if Pop.dim_in(j)>1
                        x_elem = varname2display(x_elem,a);         % copy a, after a comma if x_elem has a subscript
                    end
                    if ~isempty(out_names)
                        x_elem = [x_elem,'(',strjoin(out_names,','),')'];
                    end
                    terms{end+1} = x_elem;                  signs(end+1) = 1;       %#ok<AGROW>
                elseif ~isempty(blk_disp.rows)
                    % Matrix-valued block: the terms of element a of the output
                    if ~isequal(blk_disp.rows(a).terms,{'0'})
                        terms = [terms,blk_disp.rows(a).terms];     signs = [signs,blk_disp.rows(a).signs];     %#ok<AGROW>
                    end
                elseif blk_disp.ok && all(Pop.C{i,j}.dims==1)
                    terms = [terms,blk_disp.terms];         signs = [signs,blk_disp.signs];                     %#ok<AGROW>
                else
                    is_explicit = false;
                end
            end
            if numel(terms)>max_terms || any(cellfun(@numel,terms)>max_chars)
                is_explicit = false;
            end
            % Left-hand side (NAME x)_i,a(s) = ..., with the indices i and a if
            % there are multiple block rows and output elements
            lhs_str = [name,' ',disp_struct.xsym];
            idx_str = '';
            if p>1
                idx_str = sub_str(i);
            end
            if Pop.dim_out(i)>1
                idx_str = [idx_str,repmat(',',1,p>1),sub_str(a)];
            end
            if ~isempty(idx_str) || ~isempty(out_names)
                lhs_str = ['(',lhs_str,')',idx_str];
            end
            if ~isempty(out_names)
                lhs_str = [lhs_str,'(',strjoin(out_names,','),')'];
            end
            rows(end+1) = struct('lhs',[lhs_str,' ='],'terms',{terms},'signs',signs);       %#ok<AGROW>
        end
    end
    if is_explicit
        disp_struct.rows = rows;
        return
    end
end

% Define the blocks as the actual operators, if all of them are small enough
if ~expand || nnz(~cellfun(@isempty,Pop.C))>max_blocks
    return
end
defs = struct('lhs',{},'terms',{},'signs',{});
for i=1:p
    for j=1:d
        if isempty(Pop.C{i,j})
            continue
        end
        blk_disp = get_display(Pop.C{i,j},disp_struct.xnames{j});
        if ~blk_disp.ok
            return
        end
        lhs = ['(',blk{i,j},' ',disp_struct.xnames{j},')'];
        if ~isempty(Pop.C{i,j}.vars.out)
            out_names = cellfun(@varname2display,reshape(Pop.C{i,j}.vars.out,1,[]),'UniformOutput',false);
            lhs = [lhs,'(',strjoin(out_names,','),')'];
        else
            lhs = [blk{i,j},' ',disp_struct.xnames{j}];
        end
        defs(end+1) = struct('lhs',[lhs,' ='],'terms',{blk_disp.terms},'signs',blk_disp.signs);   %#ok<AGROW>
    end
end
disp_struct.defs = defs;

end
