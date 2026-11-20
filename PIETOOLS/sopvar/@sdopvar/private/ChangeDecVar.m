function P = ChangeDecVar(P,Zd,loc)
% P = CHANGEDECVAR(P,ZD) rewrites the 'sdopvar' P over the decision variable
% list ZD: the operator is unchanged, the rows of every params.B follow the
% order of ZD, with zero rows for variables P does not use.
%
% INPUTS
% - P:   'sdopvar' object;
% - Zd:  list of decision variables, containing every entry of P.Zd;
% - loc: (optional) vector with the position in Zd of each entry of P.Zd,
%        which skips the O(q) search of Zd. It is checked by size, range
%        and three names only, so a wrong permutation of the right length
%        is NOT caught: LOC must come from the list itself;
%
% OUTPUTS
% - P:   'sdopvar' object representing the same operator, with P.Zd = Zd,
%        the rows of each B{i} rearranged to follow Zd, and zero rows for
%        the decision variables newly introduced;
%
% NOTES
% The rows move by one product per cell, B_new = Pi*B with Pi(loc(j),j) = 1
% (Pi = T1' for d_old = T1*d_new). A name repeated in P.Zd sums its rows.
% Memory: Pi, 24 B per old decision variable.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding DB - 08/10/2026

% The lists as rows, for comparison and indexing
Zd_old = P.Zd(:).';
Zd_new = Zd(:).';
if nargin<3 || isempty(loc)
    % Same list in the same order: nothing to move
    if isequal(Zd_old,Zd_new)
        P.Zd = Zd_new;
        return
    end
    % Every old decision variable must appear in the new list
    [tf,loc] = ismember(Zd_old,Zd_new);
    if any(~tf)
        error('New decision-variable list must contain all existing variables.');
    end
else
    % A supplied map is trusted for the names - that is the saving - but a
    % map built for another list is caught by its size, its range, and a
    % name check at three positions.
    loc = loc(:).';
    nl = numel(Zd_old);
    probe = unique([1,ceil(nl/2),nl]);
    if numel(loc)~=nl || any(loc<1) || any(loc>numel(Zd_new)) ||...
            ~isequal(Zd_new(loc(probe)),Zd_old(probe))
        error('ChangeDecVar:badLoc',['LOC must give, for each entry of '...
            'P.Zd, its position in Zd.']);
    end
end
% old number of rows in each B{i}.
n_old = numel(Zd_old);
% new number of rows in each B{i}.
n_new = numel(Zd_new);
params = P.params;
% Pi = T1': row loc(j) of Pi*B is row j of B, the other rows are zero.
% Exact: each output entry is 1 times one input entry. LOC is made double,
% as sparse() refuses mixed index classes.
if n_old>0,     Pi = sparse(double(loc(:)),(1:n_old).',1,n_new,n_old); end
for ii=1:numel(params.A)
    Aii = params.A{ii};
    Bii = params.B{ii};
    % An empty B stands for no decision variable rows
    if isempty(Bii)
        Bii = sparse(0,numel(Aii));
    end
    block = numel(Aii);
    % Move the old rows to their new positions; new variables get zero rows
    if n_old>0
        Bii = reshape(Bii,n_old,block);
        % Pi*full(B) would be full: B is made sparse first
        if ~issparse(Bii),  Bii = sparse(Bii);  end
        params.B{ii} = Pi*Bii;
    else
        params.B{ii} = sparse(n_new,block);
    end
end
P.params = params;
P.Zd = Zd_new;
end
