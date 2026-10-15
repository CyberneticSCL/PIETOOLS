function P_out = repmat(P_in,varargin)
% P_OUT = REPMAT(P_IN,REPS) or REPMAT(P_IN,R,C) tiles the 'polyopvar'
% object P_IN REPS(1) times along dimension 1 (rows) and REPS(2) times
% along dimension 2 (columns).
%
% INPUTS
% - P_in:       'polyopvar' object;
% - varargin:   replication factors, as for the builtin 'repmat': either
%               REPS (positive integer, or 1 x 2 array of positive
%               integers; a scalar is used for both dimensions), or the
%               two positive integers R and C;
%
% OUTPUTS
% - P_out:      'polyopvar' object representing P_in tiled REPS times;
%
% NOTES
% - Only column replication factor 1 is supported: tiling along
%   dimension 2 would produce an object with matdim(2)>1, and
%   matrix-valued (multi-column) 'polyopvar' objects are not supported
%   (see horzcat and initialize).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - repmat
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

if nargin==2
    reps = varargin{1};
elseif nargin==3
    reps = [varargin{1},varargin{2}];
else
    error("repmat expects 2 or 3 input arguments.")
end
if isscalar(reps)
    reps = [reps,reps];
end
if numel(reps)~=2 || any(reps<1) || any(round(reps)~=reps)
    error("Replication factors must be specified as positive integers.")
end
if reps(2)~=1
    error("repmat of 'polyopvar' objects is currently only supported for column replication factor 1.")
end

P_out = P_in;
for i=2:reps(1)
    P_out = [P_out; P_in];
end

end
