function [C,Zd] = merge_dvar_lists(C,Zds,src)
% [C,ZD] = MERGE_DVAR_LISTS(C,ZDS,SRC) puts every block of the grid C onto
% one decision variable list ZD, given the lists ZDS of the operands that
% were concatenated into C and the grid SRC saying which operand each block
% came from. The operands are K >= 1 grids concatenated into C, e.g. two
% summands, or each 'sdopvar' block of a grid as its own operand.
%
% INPUTS
% - C:      block grid; 'sdopvar' blocks are on their own operand's list;
% - Zds:    1 x K cell of operand decision variable lists, as columns;
%           empty for an operand that was a promoted 'copvar';
%           empty also for sdopvar blocks with no decision variables;
% - src:    grid the size of C, the operand index of each block (0 = []);
%
% OUTPUTS
% - C:      the grid, every 'sdopvar' block on ZD and sharing that one
%           array, as the container invariant requires;
% - Zd:     the container's decision variable list, a column cellstr.
%
% NOTES
% When every operand that OWNS an 'sdopvar' block carries the same list -
% the normal case, operands declared by one call or carried through
% arithmetic together - it is reused and no block is remapped; each is
% only pointed at the one array, so concatenation stays O(blocks). Owners,
% not nonempty lists: sdopvar blocks on an EMPTY list must still be moved
% onto the union, where they get zero rows.
%
% Otherwise ONE 'unique(...,'stable')' over the concatenated lists gives
% both the union and, through its third output, the position in the union
% of every entry of every list: exactly the row map of each operand. Each
% block is then moved by the slice of that map of its operand, O(nnz(B)).
% An operand whose map is the identity - its list IS the union - keeps its
% blocks as they are.
%
% The map is TRUSTED to describe the rows of each block, which the container
% invariant guarantees and only direct assignment to the properties of a
% block can break. The block update checks the length, range and three
% names of the map: a map for a list of another length is caught, a
% permutation is not.
%
% Reconcile ONCE, before any block loop: done as a by-product of the loop,
% a block whose partner is absent passes through untouched and keeps the
% list of its operand, leaving the result holding blocks on two lists and
% breaking the container invariant.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - merge_dvar_lists
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding MMP, 09/25/2026

% Operands with a nonempty list
nz = ~cellfun(@isempty,Zds);
% Operands owning an sdopvar block, whatever their list; with none, the
% nonempty lists decide.
own = false(size(nz));
for ii = 1:numel(C)
    if isa(C{ii},'sdopvar') && src(ii)>0,   own(src(ii)) = true;    end
end
if ~any(own),   own = nz;   end
if ~any(own)
    Zd = cell(0,1);
    return
end
% Shortcut: every owner carries the list of the first owner
k1 = find(own,1);
Zd = Zds{k1};
same = true;
for k = find(own)
    if k~=k1 && ~isequal(Zds{k},Zd)
        same = false;
        break
    end
end
if same
    % Point every block at the one array; a block sum holds its own copy.
    for ii = 1:numel(C)
        if isa(C{ii},'sdopvar'),    C{ii}.Zd = Zd;  end
    end
    return
end

% The union, and every operand's row map, from one sort: cat = Zd(ic).
for k = 1:numel(Zds)
    if isempty(Zds{k}),     Zds{k} = cell(0,1);     end
end
[Zd,~,ic] = unique(vertcat(Zds{:}),'stable');
off = cumsum([0,cellfun(@numel,Zds)]);
% Each owner's map once, not per block, and whether it is the identity.
nZ = numel(Zd);     loc = cell(size(nz));   ident = false(size(nz));
for k = find(own)
    loc{k} = ic(off(k)+1:off(k+1));
    ident(k) = numel(loc{k})==nZ && isequal(loc{k}(:),(1:nZ)');
end
% Move each 'sdopvar' block onto the union by the row map of its operand
for ii = 1:numel(C)
    if ~isa(C{ii},'sdopvar'),   continue,   end
    k = src(ii);
    if ident(k)
        % already on the union, in its order
    elseif numel(C{ii}.Zd)==numel(loc{k})
        C{ii} = setdvars(C{ii},Zd,loc{k});
    else
        % Not on its operand's list, which the container invariant rules
        % out; fall back to the name search rather than trust the map.
        C{ii} = setdvars(C{ii},Zd);
    end
    C{ii}.Zd = Zd;          % share the one array, not setdvars' fresh row
end

end
