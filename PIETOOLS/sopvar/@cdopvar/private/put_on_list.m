function C = put_on_list(C,Zd)
% C = PUT_ON_LIST(C,Zd) puts every 'sdopvar' block in the cell C onto the
% decision variable list Zd, sharing the one array. 'sopvar' blocks and zero
% blocks are untouched.
%
% NO CALLER since 09/26/2026: @cdopvar/plus (09/25) and @cdopvar/mtimes     % MMP, 09/26/2026
% (09/26) now reconcile through 'merge_dvar_lists', which builds every      % MMP, 09/26/2026
% operand's row map from one sort instead of a name search per block. Kept  % MMP, 09/26/2026
% for reference.                                                            % MMP, 09/26/2026
% (was) Used by @cdopvar/plus and @cdopvar/mtimes, which both reconcile the % MMP, 09/26/2026 (was)
% (was) decision variables ONCE before their block loop. Doing it as a by-product % MMP, 09/26/2026 (was)
% Both reconciled the decision variables ONCE before their block loop.      % MMP, 09/26/2026
% Doing it as a by-product                                                  % MMP, 09/26/2026
% of the loop is wrong: a block whose partner is absent passes through
% untouched and would keep its operand's list, leaving the result holding
% blocks on two lists and breaking the container invariant.
%
% 'setdvars' returns immediately for a block already on Zd, so the O(nnz)
% row remap is only paid where it is real; the final assignment then shares
% the one array rather than the fresh Zd(:).' that 'setdvars' stores.
%
% Initial coding MMP, 09/17/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Moved from
%                  @mdopvar/private/ with the class.
% MMP, 09/26/2026: Header states it has no caller left (plus and mtimes use
%                  'merge_dvar_lists'). Doc only.

for ii = 1:numel(C)
    if isa(C{ii},'sdopvar')
        C{ii} = setdvars(C{ii},Zd);
        C{ii}.Zd = Zd;
    end
end
end
