function [CA,CB,Zd] = merge_dvar_pair(CA,CB,ZdA,ZdB)
% [CA,CB,ZD] = MERGE_DVAR_PAIR(CA,CB,ZDA,ZDB) puts the blocks of two
% operand grids onto one decision variable list ZD.
%
% INPUTS
% - CA, CB:   block grids of the two operands, each 'sdopvar' block on its
%             own operand's list;
% - ZdA, ZdB: the operands' decision variable lists (container Zd), either
%             orientation; empty for a promoted 'opvar';
%
% OUTPUTS
% - CA, CB:   the grids, every 'sdopvar' block on ZD. Unchanged when the
%             two lists are equal;
% - Zd:       the common list, a column cellstr: ZdA when the lists are
%             equal, else their stable union.
%
% NOTES
% ZD is the common decision basis d3 of the lists d1 = ZdA and d2 = ZdB,
% d1 = T1*d3 and d2 = T2*d3 with T1, T2 selection matrices. The two grids
% are merged as one, operand 1 = A and 2 = B, so one sort yields both row
% maps; when one list is empty (a promoted 'opvar'), the other list is
% kept and its blocks are only pointed at it, O(blocks).
%
% A list is reshaped to a column only when it is not one: ZD(:) on a
% q-length column copies its pointer array (12.3 ms at q = 1.07e6,
% measured), and a container stores its list as a column. 'isequal'
% compares sizes first, so unequal lengths are settled without comparing
% names; equal lists cost O(q) string comparisons.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - merge_dvar_pair
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
% Initial coding MMP, 09/30/2026

% Lists as columns; equal lists need no merge
if ~iscolumn(ZdA),  ZdA = ZdA(:);   end
if ~iscolumn(ZdB),  ZdB = ZdB(:);   end
Zd = ZdA;
if isequal(ZdA,ZdB),    return,     end
% Merge both grids as one, operand 1 = A and 2 = B, then split them again
nA = numel(CA);
src = [ones(1,nA), 2*ones(1,numel(CB))];
[Cab,Zd] = merge_dvar_lists([CA(:);CB(:)]',{ZdA,ZdB},src);
CA = reshape(Cab(1:nA),size(CA));
CB = reshape(Cab(nA+1:end),size(CB));

end
