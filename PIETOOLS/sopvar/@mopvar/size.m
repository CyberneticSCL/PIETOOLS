function varargout = size(P,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [M,N] = SIZE(P) returns the BLOCK dimensions of the 'mopvar' P: the number
% of output spaces M and input spaces N, i.e. size(P.C). SIZE(P,1) and
% SIZE(P,2) return them individually.
%
% This is the block grid, NOT the total component count - which differs from
% @sopvar/size, where the return is matrix dimensions. Component counts are
% P.dim_out (M x 1) and P.dim_in (N x 1); the concatenated spaces have
% sum(P.dim_out) and sum(P.dim_in) components. The grid is what a container
% is indexed by and what 'verify' iterates over, and with mixed L2 spaces a
% single component total is rarely the quantity wanted, so a caller who
% needs one should name it.
%
% See also VERIFY, MOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - size(mopvar)
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

blkdim = [size(P.C,1),size(P.C,2)];

if nargin==2
    varargout = {blkdim(dim)};
elseif nargout<=1
    varargout = {blkdim};
else
    varargout = {blkdim(1),blkdim(2)};
end

end
