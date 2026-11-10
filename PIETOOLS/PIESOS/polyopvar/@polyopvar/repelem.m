function P_out = repelem(P_in,R,varargin)
% P_OUT = REPELEM(P_IN,R) replicates each row i of 'polyopvar' object
% P_IN R(i) times (or R times, if R is scalar), analogous to the builtin
% 'repelem' applied along dimension 1.
%
% INPUTS
% - P_in:       'polyopvar' object;
% - R:          scalar, or nr x 1 array of nonnegative integers (nr =
%               number of rows of P_in), specifying how many times to
%               repeat each row;
% - varargin:   (optional) column replication factor C, which must be 1
%               (or empty);
%
% OUTPUTS
% - P_out:      'polyopvar' object with rows of P_in repeated as
%               specified;
%
% NOTES
% - Only row-wise replication is supported: repelem(P_in,R,C) with C~=1
%   is rejected, for the same reason column replication is rejected in
%   repmat.m.
% - At least one entry of R must be positive, otherwise an error is
%   thrown.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - repelem
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
% DJ, 10/15/2026: Initial coding;

if nargin==3
    if ~(isequal(varargin{1},1) || isempty(varargin{1}))
        error("repelem of 'polyopvar' objects is currently only supported for column replication factor 1.")
    end
elseif nargin~=2
    error("repelem expects 2 or 3 input arguments.")
end

nr = size(P_in,1);
if isscalar(R)
    R = R*ones(nr,1);
end
R = R(:);
if numel(R)~=nr || any(R<0) || any(round(R)~=R)
    error("Replication counts must be nonnegative integers, one per row of P_in.")
end

P_out = [];
for row_idx=1:nr
    Pi = subsref(P_in,substruct('()',{row_idx,':'}));
    for k=1:R(row_idx)
        if isempty(P_out)
            P_out = Pi;
        else
            P_out = [P_out; Pi];
        end
    end
end
if isempty(P_out)
    error("repelem with all-zero replication counts is not currently supported.")
end

end
