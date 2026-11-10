function disp_struct = get_display(Rop,name,expand)
% DISP_STRUCT = GET_DISPLAY(ROP,NAME,EXPAND) constructs the lines of the
% display of the block matrix ROP of tensor-PI operators, for use in the
% 'display' method of 'tensopmat',
%   NAME = [C_11 ... C_1d] : [X_1] --> [Y_1]
%          [ :   .    :  ]   [ : ]     [ : ]
%          [C_p1 ... C_pd]   [X_d]     [Y_p]
% with X_j the space that column j acts on, and Y_i the space that row i
% maps to, and the operator written out, or the definitions of the blocks,
% if they are small enough.
%
% INPUTS
% - Rop:    'tensopmat' object, with p x d blocks;
% - name:   (optional) char, the name of the block matrix. Defaults to
%           'Rop';
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
%       notes:  cellstr, the lines that list the blocks that are not of
%               type 'tensopvar';
%       defs:   struct array with the fields lhs, terms and signs, for
%               each of the blocks that is displayed as the actual
%               operator, in the order of the rows, with lhs the left-hand
%               side, e.g. '(C_11 x_1)(s) =', and terms and signs as for
%               the output disp_struct of 'tensopvar/get_display';
%       out_args:   p x 1 cell of cellstr, the names of the output
%               variables of each block row, with the copy as subscript if
%               there are multiple copies;
%       rows:   struct array with the fields lhs, terms and signs, for each
%               block row of Rop, with lhs the left-hand side, e.g.
%               '(Rop x)_1(s) =', and terms and signs for the sum of the
%               blocks in the row acting on x_j, see
%               'tensopvar/get_display' and 'intop/get_display'. Empty if
%               the operator is not displayed explicitly, see the notes;
%       xnames: 1 x d cellstr, the names x_j of the functions that the
%               blocks act on;
%       xsym:   char, the name x of the function that Rop acts on, from
%               which x_j follows;
%       is_id:  logical, true if Rop is an identity operator, see
%               'is_id_op', in which case 'display' shows (Rop x)_i(s) =
%               x_i(s) for each block row i. A block that is an identity
%               operator is defined as (C_ij x_j)(s) = x_j(s)
%               = x_j;
%
% NOTES
% - The blocks are named C_ij, dropping i if there is only one block row,
%   and j if there is only one block column. A block that is a scalar of
%   type 'double' is displayed as its value, and a zero block as 0.
% - The spaces follow from Rop.dims, Rop.dom and Rop.depmat. A space is
%   L_2^m[S] for m the dimension (dropped for m=1) and S the domain of the
%   functions, the product of the intervals of the variables that they
%   depend on, with the interval raised to the power c if the function
%   depends on c copies of the variable (e.g. from a tensor product).
%   With a single interval, S is that interval, e.g. L_2[0,1]. If the
%   functions depend on no variables, the space is R^m, or R for m=1. The
%   names of the variables are not displayed. The subscript 2 and the
%   powers are displayed as subscript and superscript characters.
% - An empty block matrix (p=0 or d=0) is displayed as NAME = [] between
%   the spaces, with the space R^0 for the absent block columns or rows,
%   and nothing is defined.
% - The operator is displayed explicitly, with a row for each block row, if
%   all blocks are scalar-valued, there are at most 5 block rows and
%   columns, every nonzero block is a 'tensopvar' with a single factor, an
%   'intop', or a scalar 'double' (in a row that maps to R), displayed
%   explicitly, and each row has at most 6 terms of at most 100 characters.
%   The blocks are then not defined.
% - The definitions of the blocks are only given if all of the blocks of
%   type 'tensopvar' with a single factor, or of type 'intop', can be
%   displayed as the actual operators, acting on x_j, with 'get_display' of
%   'tensopvar' and 'intop', and there are at most 6 of them. The
%   other blocks, of type 'tensopvar' with more than one factor, or of
%   other types, are not given. Then none of the definitions is given.

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
    name = 'Rop';
end
if nargin<3
    expand = true;
end
[p,d] = size(Rop.ops);

[~,sub_dig,sup_dig] = varname2display('0123456789');        % subscript and superscript digits
sub_str = @(n) sub_dig(sprintf('%d',n)-'0'+1);
sup_str = @(n) sup_dig(sprintf('%d',n)-'0'+1);
L2_sym = ['L',sub_str(2)];                              % L_2

blk = cell(p,d);
types = cell(p,d);
for i=1:p
    for j=1:d
        op_ij = Rop.ops{i,j};
        types{i,j} = class(op_ij);
        if isempty(op_ij)
            blk{i,j} = '0';
            continue
        elseif isa(op_ij,'double') && isscalar(op_ij)
            blk{i,j} = sprintf('%g',op_ij);
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
% rows map to (s=2), from the dimensions, domains and degrees in the variables
if ~Rop.chkval
    Rop = initialize(Rop);
end
dims = {Rop.dims.in,Rop.dims.out};
deps = {Rop.depmat.in,Rop.depmat.out};
doms = {Rop.dom.in,Rop.dom.out};
spc = cell(1,2);
out_args = cell(numel(dims{2}),1);          % names of the output variables of each row
for s=1:2
    spc{s} = cell(numel(dims{s}),1);
    for l=1:numel(dims{s})
        m = dims{s}(l);
        var_idcs = find(deps{s}(l,:)>0);
        if s==2
            % Variables of the row, with the copy as subscript if there are multiple
            out_args{l} = cell(1,0);
            for q=var_idcs
                for c=1:deps{2}(l,q)
                    if deps{2}(l,q)>1
                        out_args{l}{end+1} = varname2display(Rop.vars.out{q},c);
                    else
                        out_args{l}{end+1} = varname2display(Rop.vars.out{q});
                    end
                end
            end
        end
        if isempty(var_idcs)
            % Functions of no variables: constants
            spc_str = 'R';
            if m~=1
                spc_str = ['R',sup_str(m)];
            end
        else
            % Product of the intervals, each to the power of its number of copies
            dom_strs = cell(1,numel(var_idcs));
            for q=1:numel(var_idcs)
                dom_strs{q} = sprintf('[%g,%g]',doms{s}(var_idcs(q),1),doms{s}(var_idcs(q),2));
                if deps{s}(l,var_idcs(q))>1
                    dom_strs{q} = [dom_strs{q},sup_str(deps{s}(l,var_idcs(q)))];
                end
            end
            if isscalar(var_idcs) && deps{s}(l,var_idcs)==1
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

disp_struct = struct('lines',{lines},'notes',{cell(1,0)}, ...
    'defs',struct('lhs',{},'terms',{},'signs',{}),'xnames',{cell(1,d)},'xsym','','is_id',false,'out_args',{out_args}, ...
    'rows',struct('lhs',{},'terms',{},'signs',{}),'spc_in',{spc{1}},'spc_out',{spc{2}});

% List the blocks that are not tensor-PI operators
notes = cell(1,0);
other_types = unique(types(~cellfun(@isempty,Rop.ops) & ~strcmp(types,'tensopvar')));
for l=1:numel(other_types)
    is_type = strcmp(types,other_types{l}) & ~cellfun(@isempty,Rop.ops) & ~cellfun(@(a) isa(a,'double') && isscalar(a),Rop.ops);
    [col_idcs,row_idcs] = find(is_type.');             % in the order of the rows
    blk_names_l = arrayfun(@(a,b) blk{a,b},row_idcs,col_idcs,'UniformOutput',false);
    if isscalar(blk_names_l)
        notes{end+1} = sprintf('  The block %s is of type ''%s''.',blk_names_l{1},other_types{l});      %#ok<AGROW>
    elseif ~isempty(blk_names_l)
        notes{end+1} = sprintf('  The blocks %s are of type ''%s''.',strjoin(blk_names_l',', '),other_types{l});  %#ok<AGROW>
    end
end
disp_struct.notes = notes;

% Name of the functions on which the blocks act, x, unless this is also the
% name of a variable (ignoring a number at the end of the name), in which case
% the first of u, w, v that is not
var_bases = regexprep([Rop.vars.out,Rop.vars.in],'_?\d+$','');
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
disp_struct.is_id = is_id_op(Rop);
if is_empty
    % Nothing to define
    return
end

% Display the operator explicitly, with each element of the output as the sum
% of the blocks acting on the elements of the input, if the blocks are
% scalar-valued, there are few elements, and each element has few terms of few
% characters.
if expand && all(Rop.dims.in==1) && all(Rop.dims.out==1) && p<=max_elems && d<=max_elems
    rows = struct('lhs',{},'terms',{},'signs',{});
    is_explicit = true;
    for i=1:p
        terms = cell(1,0);      signs = zeros(1,0);
        for j=1:d
            op_ij = Rop.ops{i,j};
            if isempty(op_ij) || (isa(op_ij,'double') && ~any(op_ij(:)))
                continue
            end
            if isa(op_ij,'tensopvar') && size(op_ij.ops,2)==1
                blk_disp = get_display(op_ij,disp_struct.xnames(j));
                is_ok = blk_disp.expanded;
            elseif isa(op_ij,'intop')
                blk_disp = get_display(op_ij,disp_struct.xnames{j});
                is_ok = blk_disp.ok;
            elseif isa(op_ij,'double') && isscalar(op_ij) && isempty(out_args{i})
                % A constant times the function x_j
                const_str = [repmat([sprintf('%.4g',abs(op_ij)),' '],1,abs(op_ij)~=1),disp_struct.xnames{j}];
                blk_disp = struct('terms',{{const_str}},'signs',sign(op_ij));
                is_ok = true;
            else
                is_ok = false;
            end
            if is_ok
                terms = [terms,blk_disp.terms];     signs = [signs,blk_disp.signs];     %#ok<AGROW>
            else
                is_explicit = false;
            end
        end
        if numel(terms)>max_terms || any(cellfun(@numel,terms)>max_chars)
            is_explicit = false;
        end
        % Left-hand side (NAME x)_i(s) = ..., with the index i if there are
        % multiple block rows
        lhs_str = [name,' ',disp_struct.xsym];
        if p>1 || ~isempty(out_args{i})
            lhs_str = ['(',lhs_str,')'];
        end
        if p>1
            lhs_str = [lhs_str,sub_str(i)];
        end
        if ~isempty(out_args{i})
            lhs_str = [lhs_str,'(',strjoin(out_args{i},','),')'];
        end
        rows(end+1) = struct('lhs',[lhs_str,' ='],'terms',{terms},'signs',signs);       %#ok<AGROW>
    end
    if is_explicit
        disp_struct.rows = rows;
        return
    end
end
% Define the blocks as the actual operators, if all of them are small enough
if ~expand || nnz(~cellfun(@isempty,Rop.ops) & ~cellfun(@(a) isa(a,'double'),Rop.ops))>max_blocks
    return
end
defs = struct('lhs',{},'terms',{},'signs',{});
for i=1:p
    for j=1:d
        op_ij = Rop.ops{i,j};
        if isempty(op_ij) || isa(op_ij,'double')
            continue
        end
        if isa(op_ij,'tensopvar') && size(op_ij.ops,2)==1
            blk_disp = get_display(op_ij,disp_struct.xnames(j));
            is_ok = blk_disp.expanded;
            lhs = ['(',blk{i,j},' ',disp_struct.xnames{j},')'];
            if is_ok && ~isempty(blk_disp.out_args)
                lhs = [lhs,'(',strjoin(blk_disp.out_args,','),')'];
            end
        elseif isa(op_ij,'intop')
            blk_disp = get_display(op_ij,disp_struct.xnames{j});
            is_ok = blk_disp.ok;
            lhs = [blk{i,j},' ',disp_struct.xnames{j}];
        else
            is_ok = false;
        end
        if ~is_ok
            return
        end
        defs(end+1) = struct('lhs',[lhs,' ='],'terms',{blk_disp.terms},'signs',blk_disp.signs);   %#ok<AGROW>
    end
end
disp_struct.defs = defs;

end
