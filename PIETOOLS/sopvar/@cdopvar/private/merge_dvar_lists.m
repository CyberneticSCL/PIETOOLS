function [C,Zd] = merge_dvar_lists(C,Zds,src)
% [C,ZD] = MERGE_DVAR_LISTS(C,ZDS,SRC) puts every block of the grid C onto
% one decision variable list ZD, given the lists ZDS of the operands that
% were concatenated into C and the grid SRC saying which operand each block
% came from. This is the rule '@cdopvar/plus' applies to two summands,
% extended to K operands for horzcat, vertcat and blkdiag.
% '@cdopvar/plus' calls it too, as a two-operand grid.                      % MMP, 09/25/2026
%
% INPUTS
% - C:      block grid; 'sdopvar' blocks are on their own operand's list;
% - Zds:    1 x K cell of operand decision variable lists, as columns;
%           empty for an operand that was a promoted 'copvar';
%           empty also for sdopvar blocks with no decision variables;       % MMP, 09/25/2026
% - src:    grid the size of C, the operand index of each block (0 = []);
%
% OUTPUTS
% - C:      the grid, every 'sdopvar' block on ZD and sharing that one      % MMP, 09/25/2026
%           array, as the container invariant requires;                     % MMP, 09/25/2026
% (was) - C:      the grid, every 'sdopvar' block on ZD;                    % MMP, 09/25/2026 (was)
% - Zd:     the container's decision variable list, a column cellstr.
%
% NOTES
% When every operand that OWNS an 'sdopvar' block carries the same list -   % MMP, 09/25/2026
% the normal case, operands declared by one call or carried through         % MMP, 09/25/2026
% arithmetic together - it is reused and no block is remapped; each is      % MMP, 09/25/2026
% only pointed at the one array, so concatenation stays O(blocks). Owners,  % MMP, 09/25/2026
% not nonempty lists: sdopvar blocks on an EMPTY list must still be moved   % MMP, 09/25/2026
% onto the union, where they get zero rows.                                 % MMP, 09/25/2026
% (was) When every nonempty list is the same one - the normal case, operands % MMP, 09/25/2026 (was)
% (was) declared by one call or carried through arithmetic together - it is % MMP, 09/25/2026 (was)
% (was) reused and no block is touched, so concatenation stays O(blocks).   % MMP, 09/25/2026 (was)
%
% Otherwise ONE 'unique(...,'stable')' over the concatenated lists gives
% both the union and, through its third output, the position in the union
% of every entry of every list: exactly each operand's row map. Each block
% is then moved by its operand's slice of that map, O(nnz(B)). The earlier
% route, 'put_on_list', lets every block search the union by name; measured
% on a two-list [D, E] at q = 4e5 that per-block 'ismember' was 4.4 of the
% 5.3 s. An operand whose map is the identity - its list IS the union -     % MMP, 09/25/2026
% keeps its blocks as they are.                                             % MMP, 09/25/2026
% (was) 5.3 s. 'setdvars' spot-checks each supplied map, so a block that is not % MMP, 09/25/2026 (was)
% (was) on its operand's list is caught rather than silently misplaced.     % MMP, 09/25/2026 (was)
%
% The map is TRUSTED to describe each block's rows, which the container     % MMP, 09/25/2026
% invariant guarantees and only direct assignment to a block's properties   % MMP, 09/25/2026
% can break. 'setdvars' checks its length, range and three names: a map     % MMP, 09/25/2026
% for a list of another length is caught, a permutation is not.             % MMP, 09/25/2026
%
% Initial coding MMP, 09/25/2026
% MMP, 09/25/2026: Three changes after review, each reproduced first.
%                  (1) Sharing is decided over the operands that OWN an
%                  'sdopvar' block, not over the nonempty lists: [D, E] with
%                  E's blocks sdopvar on an EMPTY list kept D's list and left
%                  E's blocks with 0 decision rows in a container of 10.
%                  (2) An operand whose map is the identity is not rebuilt:
%                  [R2; Q1] with R2 already on the union ran 8 'setdvars'
%                  where 4 are needed. (3) On the shared-list shortcut every
%                  block is pointed at the one Zd array; a block sum carries
%                  a private copy of the list per block. The note that
%                  'setdvars' catches a misplaced block was wrong (a
%                  permutation passes) and is replaced. Now also called by
%                  '@cdopvar/plus'.

nz = ~cellfun(@isempty,Zds);
% Operands owning an sdopvar block, whatever their list; with none, the     % MMP, 09/25/2026
% nonempty lists decide, as before.                                         % MMP, 09/25/2026
own = false(size(nz));                                                      % MMP, 09/25/2026
for ii = 1:numel(C)                                                         % MMP, 09/25/2026
    if isa(C{ii},'sdopvar') && src(ii)>0,   own(src(ii)) = true;    end     % MMP, 09/25/2026
end                                                                         % MMP, 09/25/2026
if ~any(own),   own = nz;   end                                             % MMP, 09/25/2026
% if ~any(nz)                                                               % MMP, 09/25/2026 (was)
if ~any(own)                                                                % MMP, 09/25/2026
    Zd = cell(0,1);
    return
end
% k1 = find(nz,1);                                                          % MMP, 09/25/2026 (was)
k1 = find(own,1);                                                           % MMP, 09/25/2026
Zd = Zds{k1};
same = true;
% for k = find(nz)                                                          % MMP, 09/25/2026 (was)
for k = find(own)                                                           % MMP, 09/25/2026
    if k~=k1 && ~isequal(Zds{k},Zd)
        same = false;
        break
    end
end
if same
    % Point every block at the one array; a block sum holds its own copy.   % MMP, 09/25/2026
    for ii = 1:numel(C)                                                     % MMP, 09/25/2026
        if isa(C{ii},'sdopvar'),    C{ii}.Zd = Zd;  end                     % MMP, 09/25/2026
    end                                                                     % MMP, 09/25/2026
    return
end

% The union, and every operand's row map, from one sort: cat = Zd(ic).
for k = 1:numel(Zds)
    if isempty(Zds{k}),     Zds{k} = cell(0,1);     end
end
[Zd,~,ic] = unique(vertcat(Zds{:}),'stable');
off = cumsum([0,cellfun(@numel,Zds)]);
% Each owner's map once, not per block, and whether it is the identity.     % MMP, 09/25/2026
nZ = numel(Zd);     loc = cell(size(nz));   ident = false(size(nz));        % MMP, 09/25/2026
for k = find(own)                                                           % MMP, 09/25/2026
    loc{k} = ic(off(k)+1:off(k+1));                                         % MMP, 09/25/2026
    ident(k) = numel(loc{k})==nZ && isequal(loc{k}(:),(1:nZ)');             % MMP, 09/25/2026
end                                                                         % MMP, 09/25/2026
for ii = 1:numel(C)
    if ~isa(C{ii},'sdopvar'),   continue,   end
    k = src(ii);
%   loc = ic(off(k)+1:off(k+1));                                            % MMP, 09/25/2026 (was)
%   if numel(C{ii}.Zd)==numel(loc)                                          % MMP, 09/25/2026 (was)
%       C{ii} = setdvars(C{ii},Zd,loc);                                     % MMP, 09/25/2026 (was)
    if ident(k)                                                             % MMP, 09/25/2026
        % already on the union, in its order                                % MMP, 09/25/2026
    elseif numel(C{ii}.Zd)==numel(loc{k})                                   % MMP, 09/25/2026
        C{ii} = setdvars(C{ii},Zd,loc{k});                                  % MMP, 09/25/2026
    else
        % Not on its operand's list, which the container invariant rules
        % out; fall back to the name search rather than trust the map.
        C{ii} = setdvars(C{ii},Zd);
    end
    C{ii}.Zd = Zd;          % share the one array, as 'put_on_list' does
end

end
