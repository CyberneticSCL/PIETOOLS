function Pout = repmat(P,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pout = repmat(P,p,q) tiles the sdopvar object P into a p-by-q block
% pattern, so that Pout is the operator whose (i,j) block is P.
%
% INPUT
% P:        sdopvar class object
% p,q:      number of copies down and across; may also be given as [p,q],
%           or as a single scalar p meaning [p,p]
%
% OUTPUT
% Pout:     sdopvar object of dimension [p*P.dims(1), q*P.dims(2)]
%
% NOTES:
% Tiling a matrix is the same as indexing it with repeated row and column
% indices, so this is one call to @sdopvar/subsref. Nothing else is needed:
% every copy shares the monomial bases and decision variables of P.
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

% Accept repmat(P,n), repmat(P,[p,q]) and repmat(P,p,q). An operator is a
% matrix of blocks, so a third replication factor is rejected rather than
% ignored, and the factors must be scalars.
if numel(varargin)>2
    error('repmat accepts at most two dimension parameters -- see header.')
elseif isscalar(varargin)
    reps = varargin{1};
    if isscalar(reps),  reps = [reps,reps];  end
else
    reps = [varargin{1},varargin{2}];
end
if numel(reps)~=2 || ~isnumeric(reps)
    error('repmat expects the number of copies as n, [p,q] or p,q.')
end

% Indexing with repeated indices tiles the operator. Called through
% 'subsref' explicitly, since paren indexing inside a class method is not
% guaranteed to dispatch to the overload.
rowidx = repmat(1:P.dims(1),1,reps(1));
colidx = repmat(1:P.dims(2),1,reps(2));
Pout = subsref(P,substruct('()',{rowidx,colidx}));

end
