function P = setdvars(P,Zd,loc)                                             % MMP, 09/25/2026
% function P = setdvars(P,Zd)                                               % MMP, 09/25/2026 (was)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = SETDVARS(P,Zd) rewrites the 'sdopvar' P over the decision variable
% list Zd, which must contain every variable P already uses. The operator is
% unchanged: only the row indexing of params.B moves, with zero rows for
% variables P does not use.
%
% P = SETDVARS(P,Zd,LOC) does the same with the row map given: LOC(i) is    % MMP, 09/25/2026
% the position in Zd of P.Zd{i}. It skips the name search, which dominates  % MMP, 09/25/2026
% when many blocks share one list; see 'ChangeDecVar'.                      % MMP, 09/25/2026
%
% INPUTS
% - P:      an 'sdopvar' object;
% - Zd:     cellstr of decision variable names, a superset of P.Zd;
% - loc:    (optional) numel(P.Zd) positions of P.Zd in Zd;                 % MMP, 09/25/2026
%
% OUTPUTS
% - P:      the same operator with P.Zd = Zd;
%
% NOTES
% Public entry point to the private 'ChangeDecVar', needed because the
% 'copvar' container holds one decision variable list for all of its blocks
% and must be able to put a block onto it. Reimplementing the row remap in
% @copvar would duplicate logic on the decision variable axis, which is the
% axis that scales to millions.
%
% Cost: O(nnz(B)) per parameter, plus the O(q) 'ismember' over the names
% that 'ChangeDecVar' performs to locate the old rows. Returns immediately
% when P.Zd already equals Zd.
%
% See also CHANGEDECVAR, SYNC_BASIS, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - setdvars
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
% Initial coding MMP, 09/17/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change.
% MMP, 09/25/2026: Optional third input LOC, passed through to
%                  'ChangeDecVar', so the container concatenations can move a
%                  block onto a merged list without a per-block name search.
%                  Without LOC the behaviour is unchanged.

if isstring(Zd)
    Zd = cellstr(Zd);
end
% P = ChangeDecVar(P,Zd(:));                                                % MMP, 09/25/2026 (was)
if nargin<3,    loc = [];   end                                             % MMP, 09/25/2026
P = ChangeDecVar(P,Zd(:),loc);                                              % MMP, 09/25/2026

end
