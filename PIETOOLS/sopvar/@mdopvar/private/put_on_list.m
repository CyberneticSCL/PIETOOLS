function C = put_on_list(C,Zd)
% C = PUT_ON_LIST(C,Zd) puts every 'sdopvar' block in the cell C onto the
% decision variable list Zd, sharing the one array. 'sopvar' blocks and zero
% blocks are untouched.
%
% Used by @mdopvar/plus and @mdopvar/mtimes, which both reconcile the
% decision variables ONCE before their block loop. Doing it as a by-product
% of the loop is wrong: a block whose partner is absent passes through
% untouched and would keep its operand's list, leaving the result holding
% blocks on two lists and breaking the container invariant.
%
% 'setdvars' returns immediately for a block already on Zd, so the O(nnz)
% row remap is only paid where it is real; the final assignment then shares
% the one array rather than the fresh Zd(:).' that 'setdvars' stores.
%
% Initial coding MMP, 09/17/2026

for ii = 1:numel(C)
    if isa(C{ii},'sdopvar')
        C{ii} = setdvars(C{ii},Zd);
        C{ii}.Zd = Zd;
    end
end
end
