function P = ChangeDecVar(P,Zd,loc)                                         % MMP, 09/25/2026
% function P = ChangeDecVar(P,Zd)                                           % MMP, 09/25/2026 (was)
% P = CHANGEDECVAR(P,ZD) rewrites the 'sdopvar' P over the decision variable
% list ZD: the operator is unchanged, the rows of every params.B follow the
% order of ZD, with zero rows for variables P does not use.
%
% INPUTS
% - P:   'sdopvar' object;
% - Zd:  list of decision variables, containing every entry of P.Zd;
% - loc: (optional) vector with the position in Zd of each entry of P.Zd;
%        see the 09/25/2026 entry below;
%
% OUTPUTS
% - P:   'sdopvar' object representing the same operator, with P.Zd = Zd,
%        the rows of each B{i} rearranged to follow Zd, and zero rows for
%        the decision variables newly introduced;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ChangeDecVar
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
% MMP, 09/25/2026: Optional third input LOC, the position in Zd of each
%                  entry of P.Zd, for a caller that moves many blocks onto
%                  one list and has the map already (the container
%                  concatenations, from one 'unique'); it skips the O(q)
%                  'ismember' per block. LOC is checked by size, range and
%                  three names, so a wrong permutation of the right length
%                  is NOT caught: LOC must come from the list itself.
% MMP, 09/29/2026: The row remap is one product per cell, B_new = Pi*B with
%                  Pi(loc(j),j) = 1 (Pi = T1' for d_old = T1*d_new, spec
%                  sec. 8.1.4), replacing a scatter into a preallocated
%                  sparse (3-D heat build: 0.9 -> 0.25 s per large call).
%                  A name repeated in P.Zd now sums its rows, which is C(d)
%                  read by name, where the scatter kept the last row. B is
%                  made sparse first (Pi*full(B) is full) and LOC double
%                  (sparse() refuses mixed index classes). Memory: Pi, 24 B
%                  per old decision variable.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): help block in the
%                  house format, log entries condensed, body comments
%                  tightened and an example line removed; code unchanged.
%
% The lists as rows, for comparison and indexing
Zd_old = P.Zd(:).';
Zd_new = Zd(:).';
if nargin<3 || isempty(loc)                                                 % MMP, 09/25/2026
    % Same list in the same order: nothing to move                          % MMP, 09/29/2026
    if isequal(Zd_old,Zd_new)                                               % MMP, 09/29/2026
        P.Zd = Zd_new;                                                      % MMP, 09/29/2026
        return                                                              % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
    % Every old decision variable must appear in the new list               % MMP, 09/29/2026
    [tf,loc] = ismember(Zd_old,Zd_new);                                     % MMP, 09/29/2026
    if any(~tf)                                                             % MMP, 09/29/2026
        error('New decision-variable list must contain all existing variables.'); % MMP, 09/29/2026
    end                                                                     % MMP, 09/29/2026
else                                                                        % MMP, 09/25/2026
    % A supplied map is trusted for the names - that is the saving - but a  % MMP, 09/25/2026
    % map built for another list is caught by its size, its range, and a    % MMP, 09/25/2026
    % name check at three positions.                                        % MMP, 09/25/2026
    loc = loc(:).';                                                         % MMP, 09/25/2026
    nl = numel(Zd_old);                                                     % MMP, 09/25/2026
    probe = unique([1,ceil(nl/2),nl]);                                      % MMP, 09/25/2026
    if numel(loc)~=nl || any(loc<1) || any(loc>numel(Zd_new)) ||...
            ~isequal(Zd_new(loc(probe)),Zd_old(probe))                      % MMP, 09/25/2026
        error('ChangeDecVar:badLoc',['LOC must give, for each entry of '...
            'P.Zd, its position in Zd.']);                                  % MMP, 09/25/2026
    end                                                                     % MMP, 09/25/2026
end                                                                         % MMP, 09/25/2026
% old number of rows in each B{i}.
n_old = numel(Zd_old);
% new number of rows in each B{i}.
n_new = numel(Zd_new);
params = P.params;
% Pi = T1' (spec sec. 8.1.4): row loc(j) of Pi*B is row j of B, the other   % MMP, 09/29/2026
% rows are zero. Exact: each output entry is 1 times one input entry.       % MMP, 09/29/2026
if n_old>0,     Pi = sparse(double(loc(:)),(1:n_old).',1,n_new,n_old); end  % MMP, 09/29/2026
for ii=1:numel(params.A)
    Aii = params.A{ii};
    Bii = params.B{ii};
    % An empty B stands for no decision variable rows
    if isempty(Bii)
        Bii = sparse(0,numel(Aii));
    end
    block = numel(Aii);
%   % Creates a new zero sparse B matrix with one row per new decision variable. % MMP, 09/29/2026 (was)
%   Bnew = sparse(n_new,block);                                             % MMP, 09/29/2026 (was)
    % Move the old rows to their new positions; new variables get zero rows
    if n_old>0
%       Bnew(loc,:) = reshape(Bii,n_old,block);                             % MMP, 09/29/2026 (was)
        Bii = reshape(Bii,n_old,block);                                     % MMP, 09/29/2026
        % Pi*full(B) would be full; the scatter it replaces stayed sparse.  % MMP, 09/29/2026
        if ~issparse(Bii),  Bii = sparse(Bii);  end                         % MMP, 09/29/2026
        params.B{ii} = Pi*Bii;                                              % MMP, 09/29/2026
    else                                                                    % MMP, 09/29/2026
        params.B{ii} = sparse(n_new,block);                                 % MMP, 09/29/2026
    end
%   params.B{ii} = Bnew;                                                    % MMP, 09/29/2026 (was)
end
P.params = params;
P.Zd = Zd_new;
end
