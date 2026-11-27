function P = setdvars(P,Zd,loc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = SETDVARS(P,Zd) rewrites the 'sdopvar' P over the decision variable
% list Zd, which must contain every variable P already uses. The operator is
% unchanged: only the row indexing of params.B moves, with zero rows for
% variables P does not use.
%
% P = SETDVARS(P,Zd,LOC) does the same with the row map given: LOC(i) is
% the position in Zd of P.Zd{i}. It skips the name search, which dominates
% when many blocks share one list. LOC is checked by size, range and three
% names only, so it must come from the list itself.
%
% INPUTS
% - P:      an 'sdopvar' object;
% - Zd:     cellstr of decision variable names, a superset of P.Zd;
% - loc:    (optional) numel(P.Zd) positions of P.Zd in Zd;
%
% OUTPUTS
% - P:      the same operator with P.Zd = Zd;
%
% NOTES
% A 'dopvar' container holds one decision variable list for all of its
% blocks; this puts a block onto such a list.
%
% Cost: O(nnz(B)) per parameter, plus, without LOC, an O(q) search of the
% names to locate the old rows. Returns immediately when P.Zd already
% equals Zd.
%
% See also CHANGEDECVAR, SYNC_BASIS, OPVAR.
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

% Names as a column cellstr; move the rows of params.B onto Zd
if isstring(Zd)
    Zd = cellstr(Zd);
end
if nargin<3,    loc = [];   end
P = ChangeDecVar(P,Zd(:),loc);

end
