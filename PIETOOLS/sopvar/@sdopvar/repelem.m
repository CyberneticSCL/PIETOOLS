function Pout = repelem(P,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pout = repelem(P,p,q) repeats each block of the sdopvar object P, p times
% down and q times across, keeping the copies of a block adjacent.
%
% INPUT
% P:        sdopvar class object
% p,q:      number of copies of each row block and each column block; may
%           also be given as [p,q]
%
% OUTPUT
% Pout:     sdopvar object of dimension [p*P.dims(1), q*P.dims(2)]
%
% NOTES:
% Repeating the elements of a matrix is the same as indexing it with
% consecutively repeated row and column indices, so this is one call to
% @sdopvar/subsref. It differs from 'repmat' only in the order of those
% indices: repelem gives 1,1,2,2,... where repmat gives 1,2,...,1,2,...
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/06/2026

% Accept repelem(P,[p,q]) and repelem(P,p,q). The two counts are kept
% separate rather than flattened, so a per-element vector of counts reaches
% the builtin 'repelem' below and is validated there; a third count is
% rejected rather than ignored.
if numel(varargin)>2
    error('An sdopvar is two-dimensional; repelem accepts at most two replication factors.')
elseif isscalar(varargin)
    reps = varargin{1};
    if numel(reps)~=2
        error('repelem expects the counts as [p,q] when given a single argument.')
    end
    rowreps = reps(1);          colreps = reps(2);
else
    rowreps = varargin{1};      colreps = varargin{2};
end

% Indexing with consecutively repeated indices repeats the blocks. Called
% through 'subsref' explicitly, since paren indexing inside a class method
% is not guaranteed to dispatch to the overload.
rowidx = repelem(1:P.dims(1),rowreps);
colidx = repelem(1:P.dims(2),colreps);
Pout = subsref(P,substruct('()',{rowidx,colidx}));

end
