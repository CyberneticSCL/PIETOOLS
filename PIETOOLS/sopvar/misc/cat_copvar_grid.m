function [C,meta,Zds,src] = cat_copvar_grid(dir,args,cls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,META,ZDS,SRC] = CAT_COPVAR_GRID(DIR,ARGS,CLS) assembles the block grid and
% metadata of a horizontal, vertical or block-diagonal concatenation of
% 'copvar' or 'cdopvar' containers.
%
% INPUTS
% - dir:    'h' for [A, B, ...], 'v' for [A; B; ...], 'd' for
%           blkdiag(A, B, ...);
% - args:   1 x K cell of operands: containers of either class, or single
%           'sopvar'/'sdopvar' blocks, which are taken as 1 x 1 containers.
%           Numeric [] operands, as MATLAB passes for [A, []], are dropped,
%           and so are containers with an empty 0 x 0 grid;
% - cls:    class of the result, 'copvar' or 'cdopvar'. Every operand is
%           promoted to it;
%
% OUTPUTS
% - C:      the concatenated block grid; zero blocks stay [];
% - meta:   struct with fields vars, dom, space_out, space_in, dim_out,
%           dim_in: the container metadata, without a decision variable
%           list;
% - Zds:    1 x K' cell holding the decision variable list of each
%           surviving operand, as a column; the stored list itself, not a
%           copy. Empty lists for 'copvar';
% - src:    grid the size of C: which surviving operand each block came
%           from, 0 for a [] block. With Zds, every block of one operand can
%           be moved by ONE row map when the decision lists are merged.
%
% NOTES
% A container grid keeps its spaces SEPARATE: [A, B] has the columns of A
% followed by the columns of B, and each block is left as it was; nothing
% is reordered.
%
% [A, B] needs A and B to agree on every ROW: the same number of block rows,
% and in row i the same output space (as a SET of variable names) and the
% same component count. [A; B] needs the same of COLUMNS. blkdiag needs
% nothing. The operands may be built over different variable registries -
% an R^n -> R^n block has an empty one - so the registries are merged, a
% variable given two different domains is an error, and every space mask is
% remapped onto the merged registry. That is what lets [Dzw, Cz] join a
% finite-dimensional block to one over L_2[s].
%
% The grid and metadata logic is the same for both container classes. The
% decision variable lists are not merged here; ZDS and SRC carry what the
% merge needs.
%
% Cost: O(total blocks) cell moves, plus set operations on the variable
% registries, which have one entry per spatial direction. Nothing here
% touches a coefficient matrix or scales with the number of decision
% variables; promoting a single 'sdopvar' block to a container is O(q) once.
%
% See also HORZCAT, VERTCAT, BLKDIAG, DERIVE_COPVAR_META, COPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cat_copvar_grid
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
% Initial coding MMP, 09/25/2026

% Check the direction; the name of the operation, for messages
if ~any(strcmp(dir,{'h','v','d'}))
    error('cat_copvar_grid:badDir',"'dir' must be 'h', 'v' or 'd'.")
end
opname = struct('h','horzcat','v','vertcat','d','blkdiag');
opname = opname.(dir);

% % % Promote every operand to the result class, dropping empties.
keep = true(1,numel(args));
for k = 1:numel(args)
    a = args{k};
    % Numeric [] operands are dropped
    if isnumeric(a) && isempty(a)
        keep(k) = false;
        continue
    end
    % A block, or a 'copvar' for a 'cdopvar' result, is promoted; any other
    % class is refused
    if strcmp(cls,'copvar')
        if isa(a,'sopvar') && ~isa(a,'sdopvar')
            a = copvar({a});
        elseif ~isa(a,'copvar')
            error([cls ':' opname 'BadOperand'],['Cannot %s a ''copvar'' with a '...
                '''%s''. Operands must be copvar/cdopvar containers or '...
                'sopvar/sdopvar blocks; convert others first, e.g. with '...
                'mat2copvar_sop.'],opname,class(a))
        end
    else
        % 'cdopvar' result: a 'copvar' is promoted, a block wrapped
        if isa(a,'copvar')
            a = cdopvar(a);
        elseif isa(a,'sopvar') || isa(a,'sdopvar')
            a = cdopvar({a});
        elseif ~isa(a,'cdopvar')
            error([cls ':' opname 'BadOperand'],['Cannot %s a ''cdopvar'' with a '...
                '''%s''. Operands must be copvar/cdopvar containers or '...
                'sopvar/sdopvar blocks; convert others first, e.g. with '...
                'mat2copvar_sop.'],opname,class(a))
        end
    end
    % An empty 0 x 0 container is the identity of concatenation.
    if numel(a.C)==0 && size(a.C,1)==0 && size(a.C,2)==0
        keep(k) = false;
        continue
    end
    args{k} = a;
end
% No operand left: the empty container
args = args(keep);
K = numel(args);
if K==0
    C = cell(0,0);
    meta = struct('vars',{cell(1,0)},'dom',zeros(0,2),...
        'space_out',false(0,0),'space_in',false(0,0),...
        'dim_out',zeros(0,1),'dim_in',zeros(0,1));
    Zds = cell(1,0);
    src = zeros(0,0);
    return
end

% % % Merge the registries: every operand comes back on the sorted union
% % % registry, its masks remapped by name; a variable given two domains is
% % % an error (cls:domConflict). Registries are short (one entry per
% % % spatial direction), so this is on the small axis.
[args{1:K}] = merge_copvar_registry(cls,opname,args{:});
vars = args{1}.vars;        dom = args{1}.dom;
S_out = cell(1,K);      S_in = cell(1,K);
D_out = cell(1,K);      D_in = cell(1,K);
for k = 1:K
    a = args{k};
    S_out{k} = a.space_out;         S_in{k} = a.space_in;
    D_out{k} = a.dim_out(:);        D_in{k} = a.dim_in(:);
end

% % % Agreement checks, and the grid.
switch dir
    case 'h'
        M = size(args{1}.C,1);
        for k = 2:K
            if size(args{k}.C,1)~=M
                error([cls ':horzcatRowMismatch'],['Operand %d of horzcat has '...
                    '%d block rows; operand 1 has %d.'],k,size(args{k}.C,1),M)
            end
            check_side(S_out{1},S_out{k},D_out{1},D_out{k},vars,k,'row',...
                'output',cls,'horzcat');
        end
        C = cell(M,0);      src = zeros(M,0);
        for k = 1:K
            C = [C, args{k}.C];                                     %#ok<AGROW>
            src = [src, k*ones(size(args{k}.C))];                   %#ok<AGROW>
        end
        space_out = S_out{1};       dim_out = D_out{1};
        space_in = vertcat(S_in{:});    dim_in = vertcat(D_in{:});
    case 'v'
        N = size(args{1}.C,2);
        for k = 2:K
            if size(args{k}.C,2)~=N
                error([cls ':vertcatColMismatch'],['Operand %d of vertcat has '...
                    '%d block columns; operand 1 has %d.'],k,size(args{k}.C,2),N)
            end
            check_side(S_in{1},S_in{k},D_in{1},D_in{k},vars,k,'column',...
                'input',cls,'vertcat');
        end
        C = cell(0,N);      src = zeros(0,N);
        for k = 1:K
            C = [C; args{k}.C];                                     %#ok<AGROW>
            src = [src; k*ones(size(args{k}.C))];                   %#ok<AGROW>
        end
        space_in = S_in{1};         dim_in = D_in{1};
        space_out = vertcat(S_out{:});  dim_out = vertcat(D_out{:});
    case 'd'
        Ms = cellfun(@(a) size(a.C,1),args);
        Ns = cellfun(@(a) size(a.C,2),args);
        C = cell(sum(Ms),sum(Ns));      src = zeros(sum(Ms),sum(Ns));
        r0 = 0;     c0 = 0;
        for k = 1:K
            C(r0+(1:Ms(k)),c0+(1:Ns(k))) = args{k}.C;
            src(r0+(1:Ms(k)),c0+(1:Ns(k))) = k;
            r0 = r0+Ms(k);      c0 = c0+Ns(k);
        end
        space_out = vertcat(S_out{:});  dim_out = vertcat(D_out{:});
        space_in = vertcat(S_in{:});    dim_in = vertcat(D_in{:});
end

% Metadata of the result; src is 0 at [] blocks
meta = struct('vars',{vars},'dom',dom,'space_out',space_out,...
    'space_in',space_in,'dim_out',dim_out,'dim_in',dim_in);

src(cellfun(@isempty,C)) = 0;

% The stored lists themselves: a container keeps its list as a column, and
% '(:)' on a q-long cell copies its pointer array, 8 bytes a variable -
% measured 6.1 MB for two operands at q = 4e5, on the path meant to be
% flat in q. Reshape only a list that is not already a column.
Zds = cell(1,K);
for k = 1:K
    z = cell(0,1);
    if strcmp(cls,'cdopvar'),   z = args{k}.Zd;     end
    if ~iscolumn(z),            z = z(:);           end
    Zds{k} = z;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function check_side(S1,Sk,D1,Dk,vars,k,what,side,cls,opname)
% CHECK_SIDE(S1,SK,D1,DK,VARS,K,WHAT,SIDE,CLS,OPNAME) errors unless the rows
% (horzcat) or columns (vertcat) of operand K match operand 1 one for one:
% the same space as a set of variables, and the same component count. The
% masks are already over the merged registry, so mask equality is set
% equality.
% INPUTS  S1, Sk: space masks of operands 1 and k; D1, Dk: their component
%         counts; vars: the merged registry; k: operand index; what: 'row'
%         or 'column'; side: 'output' or 'input'; cls, opname: for the
%         error identifier and message.

for r = 1:size(S1,1)
    if ~isequal(S1(r,:),Sk(r,:))
        error([cls ':' opname 'SpaceMismatch'],['Operand %d of %s maps %s '...
            '{%s} in block %s %d, but operand 1 maps %s {%s}.'],k,opname,...
            ternary(strcmp(side,'output'),'into','out of'),...
            strjoin(vars(Sk(r,:)),','),what,r,...
            ternary(strcmp(side,'output'),'into','out of'),...
            strjoin(vars(S1(r,:)),','))
    end
    if D1(r)~=Dk(r)
        error([cls ':' opname 'DimMismatch'],['Operand %d of %s has %s '...
            'dimension %d in block %s %d, but operand 1 has %d.'],...
            k,opname,side,Dk(r),what,r,D1(r))
    end
end

end


function s = ternary(c,a,b)
% S = TERNARY(C,A,B) is A when C is true and B otherwise.
% INPUTS  c: logical; a, b: any values.  OUTPUT  s: a or b.
if c,   s = a;  else,   s = b;  end
end
