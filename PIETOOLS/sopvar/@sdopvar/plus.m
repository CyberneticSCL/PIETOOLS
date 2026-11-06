function C = plus(A,B)
% C = PLUS(A,B) returns the sum A + B of two 'sdopvar' objects (or of an
% 'sdopvar' and a fixed 'sopvar' object), affine in the decision variables.
%
% INPUTS
% - A, B:   'sdopvar' objects of the same size, mapping between the same
%           spaces, or one of them a fixed 'sopvar' (promoted onto the
%           decision variables of the other), or a 'dpvar' (scalar*I or
%           matrix multiplier, via 'dpvar_op_copvar'). Text is refused;
%
% OUTPUTS
% - C:      'sdopvar' object representing A + B, on the merged decision-
%           variable and monomial bases; computed by 'plus_batch';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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
% MMP, 09/07/2026: Accept a fixed 'sopvar' summand, promoted with
%                  'sopvar2sdopvar'; both A+Pop and Pop+A dispatch here.
% MMP, 09/21/2026: Variable lists compared with 'isequal', not
%                  'any(~strcmp(...))', which is false for lists of
%                  different length ('{}' against a 1 x 0 cellstr).
% MMP, 09/25/2026: Container classes renamed mopvar -> copvar, mdopvar ->
%                  cdopvar; renamed here: posmopvar -> poscopvar.
% MMP, 09/26/2026: A 'sopvar' operand promoted onto the other's nonempty
%                  list lets 'sync_basis' skip its O(q) list comparison
%                  ('shared_Zd').
% MMP, 09/29/2026: A dpvar summand is routed to 'dpvar_op_copvar' (a scalar
%                  is scalar*I on a square block, a matrix of the block's
%                  size its multiplier); before, a dpvar never dispatched
%                  here. Numeric summands are not routed, as before.
% MMP, 09/29/2026: The sum is 'plus_batch' with two operands, which runs the
%                  four stages this file repeated for N = 2 (promotion,
%                  checks, 'sync_basis', the 'apply_basis_map' sum) and had
%                  drifted from it; that code is deleted here (BEGIN/END
%                  below), so the 09/07, 09/21 and 09/26/2026 entries above
%                  describe code in 'plus_batch'. Error texts are now
%                  plus_batch's. A text summand is refused here, since
%                  plus_batch reads text as its option.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): what-line, log
%                  entries condensed to defect and fix, the deletion note in
%                  the body shortened, a (was) line of a comment removed;
%                  code unchanged.
%

% dpvar summand: legacy scalar*I / matrix multiplier, before                % MMP, 09/29/2026
% the call to 'plus_batch' below.                                           % MMP, 09/29/2026
if isa(A,'dpvar') || isa(B,'dpvar')                                         % MMP, 09/29/2026
    C = dpvar_op_copvar('plus',A,B);                                        % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
% BEGIN MMP, 09/29/2026: delegation to 'plus_batch', which holds the checks, % MMP, 09/29/2026
% the synchronization and the sum this file had for two operands, and       % MMP, 09/29/2026
% promotes a fixed 'sopvar' operand onto the list Zd_p of its first         % MMP, 09/29/2026
% 'sdopvar' operand; Zd_p below only sets the 'shared_Zd' flag.             % MMP, 09/29/2026
% Text is plus_batch's option, never a summand: refuse it here.             % MMP, 09/29/2026
if ischar(A) || isstring(A) || ischar(B) || isstring(B)                     % MMP, 09/29/2026
    error("A summand of an 'sdopvar' must be an operator, not text.")       % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
% A nonempty Zd_p is stored unchanged by 'sopvar2sdopvar', so both lists    % MMP, 09/26/2026
% are then the same array. An empty one is replaced by cell(0,1), and the   % MMP, 09/26/2026
% comparison is O(1) there anyway, so it is left to sync_basis.             % MMP, 09/26/2026
shared_Zd = false;                                                          % MMP, 09/26/2026
if isa(A,'sopvar') || isa(B,'sopvar')                                       % MMP, 09/07/2026
    if      isa(B,'sdopvar'),   Zd_p = B.Zd;                                % MMP, 09/07/2026
    elseif  isa(A,'sdopvar'),   Zd_p = A.Zd;                                % MMP, 09/07/2026
    else,                       Zd_p = cell(0,1);                           % MMP, 09/07/2026
    end                                                                     % MMP, 09/07/2026
%   if isa(A,'sopvar'),     A = sopvar2sdopvar(A,Zd_p);     end             % MMP, 09/07/2026 % MMP, 09/29/2026 (was)
%   if isa(B,'sopvar'),     B = sopvar2sdopvar(B,Zd_p);     end             % MMP, 09/07/2026 % MMP, 09/29/2026 (was)
    shared_Zd = ~isempty(Zd_p);                                             % MMP, 09/26/2026
end                                                                         % MMP, 09/07/2026
if shared_Zd                                                                % MMP, 09/29/2026
    C = plus_batch(A,B,'shared_Zd');                                        % MMP, 09/29/2026
else                                                                        % MMP, 09/29/2026
    C = plus_batch(A,B);                                                    % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
% END MMP, 09/29/2026
end
