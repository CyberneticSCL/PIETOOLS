function [C,Zd] = merge_dvar_lists(C,Zds,src)
% [C,ZD] = MERGE_DVAR_LISTS(C,ZDS,SRC) puts every block of the grid C onto
% one decision variable list ZD, given the lists ZDS of the operands that
% were concatenated into C and the grid SRC saying which operand each block
% came from. This is the rule '@cdopvar/plus' applies to two summands,
% extended to K operands for horzcat, vertcat and blkdiag.
%
% INPUTS
% - C:      block grid; 'sdopvar' blocks are on their own operand's list;
% - Zds:    1 x K cell of operand decision variable lists, as columns;
%           empty for an operand that was a promoted 'copvar';
% - src:    grid the size of C, the operand index of each block (0 = []);
%
% OUTPUTS
% - C:      the grid, every 'sdopvar' block on ZD;
% - Zd:     the container's decision variable list, a column cellstr.
%
% NOTES
% When every nonempty list is the same one - the normal case, operands
% declared by one call or carried through arithmetic together - it is
% reused and no block is touched, so concatenation stays O(blocks).
%
% Otherwise ONE 'unique(...,'stable')' over the concatenated lists gives
% both the union and, through its third output, the position in the union
% of every entry of every list: exactly each operand's row map. Each block
% is then moved by its operand's slice of that map, O(nnz(B)). The earlier
% route, 'put_on_list', lets every block search the union by name; measured
% on a two-list [D, E] at q = 4e5 that per-block 'ismember' was 4.4 of the
% 5.3 s. 'setdvars' spot-checks each supplied map, so a block that is not
% on its operand's list is caught rather than silently misplaced.
%
% Initial coding MMP, 09/25/2026

nz = ~cellfun(@isempty,Zds);
if ~any(nz)
    Zd = cell(0,1);
    return
end
k1 = find(nz,1);
Zd = Zds{k1};
same = true;
for k = find(nz)
    if k~=k1 && ~isequal(Zds{k},Zd)
        same = false;
        break
    end
end
if same
    return
end

% The union, and every operand's row map, from one sort: cat = Zd(ic).
for k = 1:numel(Zds)
    if isempty(Zds{k}),     Zds{k} = cell(0,1);     end
end
[Zd,~,ic] = unique(vertcat(Zds{:}),'stable');
off = cumsum([0,cellfun(@numel,Zds)]);
for ii = 1:numel(C)
    if ~isa(C{ii},'sdopvar'),   continue,   end
    k = src(ii);
    loc = ic(off(k)+1:off(k+1));
    if numel(C{ii}.Zd)==numel(loc)
        C{ii} = setdvars(C{ii},Zd,loc);
    else
        % Not on its operand's list, which the container invariant rules
        % out; fall back to the name search rather than trust the map.
        C{ii} = setdvars(C{ii},Zd);
    end
    C{ii}.Zd = Zd;          % share the one array, as 'put_on_list' does
end

end
