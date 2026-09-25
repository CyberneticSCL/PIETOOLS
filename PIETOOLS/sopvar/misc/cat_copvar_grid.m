function [C,meta,Zds,src] = cat_copvar_grid(dir,args,cls)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [C,META,ZDS,SRC] = CAT_COPVAR_GRID(DIR,ARGS,CLS) assembles the block grid and
% metadata of a horizontal, vertical or block-diagonal concatenation of
% containers, for the 'horzcat', 'vertcat' and 'blkdiag' methods of 'copvar'
% and 'cdopvar'.
%
% INPUTS
% - dir:    'h' for [A, B, ...], 'v' for [A; B; ...], 'd' for
%           blkdiag(A, B, ...);
% - args:   1 x K cell of operands: containers of either class, or single
%           'sopvar'/'sdopvar' blocks, which are taken as 1 x 1 containers.
%           Numeric [] operands, as MATLAB passes for [A, []], are dropped,
%           and so are containers with an empty 0 x 0 grid;
% - cls:    class of the result, 'copvar' or 'cdopvar'. Every operand is
%           promoted to it; the caller decides it;
%
% OUTPUTS
% - C:      the concatenated block grid; zero blocks stay [];
% - meta:   struct with fields vars, dom, space_out, space_in, dim_out,
%           dim_in, for the two-argument constructor. No decision variable
%           list: a 'cdopvar' caller reconciles that itself;
% - Zds:    1 x K' cell holding each surviving operand's decision variable
%           list, for that reconciliation, as a column; the stored list
%           itself, not a copy. Empty lists for 'copvar';
% - src:    grid the size of C: which surviving operand each block came
%           from, 0 for a [] block. With Zds it lets the reconciliation
%           move every block of one operand by ONE row map.
%
% NOTES
% A container grid keeps its spaces SEPARATE: [A, B] has the columns of A
% followed by the columns of B, and each block is left as it was. This is
% not what 'dopvar/horzcat' does - it merges the finite-dimensional parts of
% all operands into one R^n block and the L_2 parts into one L_2 block - and
% it is why the container needs no 'ambiguous ordering' warning: nothing is
% reordered. Positive operators built over the result take its space list
% directly (see 'poscopvar').
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
% This routine is shared by both container classes, like
% 'derive_copvar_meta', because the grid and metadata logic is identical
% for the two and is the part that would silently diverge if duplicated.
% The decision variable list is not handled here, since reconciling it
% needs @cdopvar's private 'merge_dvar_lists', which uses SRC below.
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

if ~any(strcmp(dir,{'h','v','d'}))
    error('cat_copvar_grid:badDir',"'dir' must be 'h', 'v' or 'd'.")
end
opname = struct('h','horzcat','v','vertcat','d','blkdiag');
opname = opname.(dir);

% % % Promote every operand to the result class, dropping empties.
keep = true(1,numel(args));
for k = 1:numel(args)
    a = args{k};
    if isnumeric(a) && isempty(a)
        keep(k) = false;
        continue
    end
    if strcmp(cls,'copvar')
        if isa(a,'sopvar') && ~isa(a,'sdopvar')
            a = copvar({a});
        elseif ~isa(a,'copvar')
            error([cls ':' opname 'BadOperand'],['Cannot %s a ''copvar'' with a '...
                '''%s''. Operands must be copvar/cdopvar containers or '...
                'sopvar/sdopvar blocks; convert others first, e.g. with '...
                'opvar2copvar.'],opname,class(a))
        end
    else
        if isa(a,'copvar')
            a = cdopvar(a);
        elseif isa(a,'sopvar') || isa(a,'sdopvar')
            a = cdopvar({a});
        elseif ~isa(a,'cdopvar')
            error([cls ':' opname 'BadOperand'],['Cannot %s a ''cdopvar'' with a '...
                '''%s''. Operands must be copvar/cdopvar containers or '...
                'sopvar/sdopvar blocks; convert others first, e.g. with '...
                'opvar2copvar.'],opname,class(a))
        end
    end
    % An empty 0 x 0 container is the identity of concatenation.
    if numel(a.C)==0 && size(a.C,1)==0 && size(a.C,2)==0
        keep(k) = false;
        continue
    end
    args{k} = a;
end
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

% % % Merge the variable registries. Each is sorted and short (one entry
% % % per spatial direction), so this is on the small axis.
allv = cell(1,0);
for k = 1:K
    allv = [allv, reshape(args{k}.vars,1,[])];                     %#ok<AGROW>
end
vars = unique(allv);                        % sorted, as the registry must be
vars = reshape(vars,1,[]);
dom = nan(numel(vars),2);
for k = 1:K
    [~,loc] = ismember(args{k}.vars,vars);
    for r = 1:numel(loc)
        d = args{k}.dom(r,:);
        if all(isnan(dom(loc(r),:)))
            dom(loc(r),:) = d;
        elseif any(dom(loc(r),:)~=d)
            error([cls ':domConflict'],['Variable ''%s'' is on [%g,%g] in '...
                'operand %d of %s and on [%g,%g] in an earlier operand.'],...
                vars{loc(r)},d(1),d(2),k,opname,dom(loc(r),1),dom(loc(r),2))
        end
    end
end
nv = numel(vars);

% % % Remap each operand's space masks onto the merged registry.
S_out = cell(1,K);      S_in = cell(1,K);
D_out = cell(1,K);      D_in = cell(1,K);
for k = 1:K
    a = args{k};
    [~,loc] = ismember(a.vars,vars);
    so = false(size(a.C,1),nv);     so(:,loc) = a.space_out;
    si = false(size(a.C,2),nv);     si(:,loc) = a.space_in;
    S_out{k} = so;                  S_in{k} = si;
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
% Rows (horzcat) or columns (vertcat) must match one for one: the same
% space as a set of variables, and the same component count. The masks are
% already over the merged registry, so mask equality is set equality.

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
if c,   s = a;  else,   s = b;  end
end
