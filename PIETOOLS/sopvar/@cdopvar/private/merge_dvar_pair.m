function [CA,CB,Zd] = merge_dvar_pair(CA,CB,ZdA,ZdB)
% [CA,CB,ZD] = MERGE_DVAR_PAIR(CA,CB,ZDA,ZDB) puts the blocks of two
% operand grids onto one decision variable list ZD, for '@cdopvar/plus' and
% '@cdopvar/mtimes', which reconcile their two operands once before the
% block loop (see 'merge_dvar_lists').
%
% INPUTS
% - CA, CB:   block grids of the two operands, each 'sdopvar' block on its
%             own operand's list;
% - ZdA, ZdB: the operands' decision variable lists (container Zd), either
%             orientation; empty for a promoted 'copvar';
%
% OUTPUTS
% - CA, CB:   the grids, every 'sdopvar' block on ZD. Unchanged when the
%             two lists are equal;
% - Zd:       the common list, a column cellstr: ZdA when the lists are
%             equal, else the union from 'merge_dvar_lists'.
%
% NOTES
% This is the rule of sopvar_implementation_notes.pdf Sec. 8.1.4 (Combine
% Basis): the common decision basis d3 with d1 = T1*d3, d2 = T2*d3. The two
% grids are handed to 'merge_dvar_lists' as one, operand 1 = A and
% 2 = B, so its one sort yields both row maps; with the fixed factor of
% 'mtimes' (empty promoted list) its owner shortcut keeps the decision
% factor's list and only points the blocks at it, O(blocks).
%
% A list is reshaped to a column only when it is not one: ZD(:) on a
% q-length column copies its pointer array (12.3 ms at q = 1.07e6,
% measured), and a container stores its list as a column. 'isequal'
% compares sizes first, so unequal lengths are settled without comparing
% names; equal lists cost O(q) string comparisons.
%
% Initial coding MMP, 09/30/2026. The two-operand wrapper of
%                'merge_dvar_lists' that @cdopvar/plus (09/25/2026) and
%                @cdopvar/mtimes (09/26/2026) each carried inline, 7 lines
%                apiece, with every list read as Zd(:): up to 5 q-length
%                copies per call. Same lists, same union, same blocks.

if ~iscolumn(ZdA),  ZdA = ZdA(:);   end
if ~iscolumn(ZdB),  ZdB = ZdB(:);   end
Zd = ZdA;
if isequal(ZdA,ZdB),    return,     end
nA = numel(CA);
src = [ones(1,nA), 2*ones(1,numel(CB))];
[Cab,Zd] = merge_dvar_lists([CA(:);CB(:)]',{ZdA,ZdB},src);
CA = reshape(Cab(1:nA),size(CA));
CB = reshape(Cab(nA+1:end),size(CB));

end
