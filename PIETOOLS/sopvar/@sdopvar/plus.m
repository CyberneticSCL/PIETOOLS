function C = plus(A,B)
% Add two sdopvar operators, preserving affine decision representation.
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
% MMP, 09/07/2026: Accept a fixed 'sopvar' summand, promoting it with
% 'sopvar2sdopvar'. Adding a known block to a decision operator is normal
% usage, and since 'sdopvar' now lists ?sopvar among its InferiorClasses
% both A+Pop and Pop+A are dispatched here.
% MMP, 09/21/2026: Compare the variable lists with 'isequal' rather than
% 'any(~strcmp(...))'. '@sopvar/plus' took the same fix on 09/07/2026 and
% this copy was missed. Reached by adding the strict-positivity identity to a
% 'poscopvar' container, whose R^n block has an empty variable list: '{}' and
% a 1 x 0 cellstr name the same space, and strcmp of two different-sized
% cellstr does not compare them.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
% mdopvar -> cdopvar, with every file and function named after them. Mechanical
% rename, no functional change. Renamed here: posmopvar -> poscopvar.
% MMP, 09/26/2026: After promoting a 'sopvar' operand onto the other's
% nonempty list, tell 'sync_basis' the lists are equal (its 'shared_Zd')
% instead of letting it prove so with an O(q) isequal. 7 of the 11 calls in
% the 2-D container Hinf build (q = 3.8e5 per list, ~8 ms each); with the
% cheaper comparison in 'sync_basis', plus there 0.21 -> 0.11 s. Same output.
% MMP, 09/29/2026: A dpvar summand is routed at the top to
% 'dpvar_op_copvar' (legacy @dopvar/plus semantics): a scalar is scalar*I on
% a square block, a matrix of the block's size its multiplier. Before, a
% dpvar never dispatched here. Numeric summands are not routed, as before.
% MMP, 09/29/2026: The sum is now 'plus_batch' with two operands. It runs
% the same four stages this file repeated for N = 2 (sopvar promotion, the
% compatibility checks, 'sync_basis', the 'apply_basis_map' sum; spec sec.
% 8.1.3-8.1.4), and the copies had drifted: this one lacked plus_batch's
% skip of summands without entries. Kept here: the dpvar routing and the
% 'shared_Zd' rule of 09/26/2026, passed as plus_batch's option. Programs of
% the 1-D, 2-D and 3-D benchmark builds unchanged bit for bit. Error texts
% are now plus_batch's ('Summands map between different spaces' etc.), and
% a numeric or polynomial summand fails with 'plus_batch is supported only
% for sdopvar and sopvar objects' rather than on its missing 'dims'. The
% 09/07/2026, 09/21/2026 and 09/26/2026 entries above now describe code in
% 'plus_batch'; that code is deleted here, see BEGIN/END below. Cost,
% measured: about 20 us more per call at one spatial variable and q <= 1e4
% (+3-8%), from plus_batch's operand handling; faster from two variables or
% q >= 1e5 (47x for a 'sopvar' summand at three variables, q = 3e6). A text
% summand is refused here, since plus_batch reads text as its option.
%

% dpvar summand: legacy scalar*I / matrix multiplier, before                % MMP, 09/29/2026
% the promotion and checks below.                                           % MMP, 09/29/2026 (was)
% the call to 'plus_batch' below.                                           % MMP, 09/29/2026
if isa(A,'dpvar') || isa(B,'dpvar')                                         % MMP, 09/29/2026
    C = dpvar_op_copvar('plus',A,B);                                        % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
% BEGIN MMP, 09/29/2026: delegation to 'plus_batch'. Deleted here, done     % MMP, 09/29/2026
% there: the dims, space, domain and term-count checks (space check         % MMP, 09/29/2026
% 09/21/2026), the 'sync_basis' call (09/07/2026, flag 09/26/2026), the     % MMP, 09/29/2026
% 'apply_basis_map' sum loop (09/07/2026), the '(was)' lines of the         % MMP, 09/29/2026
% pairwise CombineDecisionBasis/UnionBasisMonomials prologue (09/07/2026),  % MMP, 09/29/2026
% and the unmarked comments of the original sum (spec sec. 8.1.3) with 2    % MMP, 09/29/2026
% commented-out lines, and the 09/07/2026 promotion comment (it named the   % MMP, 09/29/2026
% uncalled 'CombineDecisionBasis'). The promotion is commented out below.   % MMP, 09/29/2026
% 'plus_batch' promotes a fixed 'sopvar' operand onto Zd_p, the list of     % MMP, 09/29/2026
% its first 'sdopvar' operand; with two operands that is the Zd_p below,    % MMP, 09/29/2026
% which only sets the flag.                                                 % MMP, 09/29/2026
% Text is plus_batch's option, never a summand: refuse it here, as the old  % MMP, 09/29/2026
% body did (it failed on the missing 'dims').                               % MMP, 09/29/2026
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
