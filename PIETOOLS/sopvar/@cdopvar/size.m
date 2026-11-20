function varargout = size(P,dim,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m,n] = SIZE(P) returns the row and column dimensions of the operator
% represented by the 'cdopvar' P: the total numbers of output and input
% components m = sum(P.dim_out) and n = sum(P.dim_in). SIZE(P,'block')
% instead returns the block grid
% [M,N] = size(P.C), the numbers of output and input spaces, and
% SIZE(P,'split') the number of components of each of those spaces.
%
% INPUTS
% - P:      'cdopvar' container;
% - dim:    (optional) 1 or 2, return only the row (1) or column (2) size;
% - opts:   (optional) char specifying what size to return. One of:
%           'tot':   Default, the total row/column dimension of the
%                    operator, sum(P.dim_out) and sum(P.dim_in);
%           'block': The number of operator blocks along the rows/columns,
%                    size(P.C,1) and size(P.C,2);
%           'split': The row/column dimension of each block row/column,
%                    P.dim_out (M x 1) and P.dim_in (N x 1);
%           'dim' and 'opts' may be passed in either order.
%
% OUTPUTS
% - varargout:  with dim, the size along that dimension only: a scalar for
%               'tot' and 'block', the M x 1 or N x 1 array for 'split'.
%               Without dim, two outputs return the row and column sizes
%               separately, and one output returns the 1 x 2 array [m,n]
%               for 'tot' and 'block', and the 1 x 2 cell
%               {P.dim_out,P.dim_in} for 'split', as the two arrays may
%               differ in length;
%
% See also VERIFY, CDOPVAR, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - size(cdopvar)
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

% Default to the total size; accept size(P,opts) and size(P,opts,dim) as
% well as size(P,dim) and size(P,dim,opts).
if nargin<3
    opts = 'tot';
end
if nargin<2
    dim = [];
end
if ischar(dim) || isstring(dim)
    if nargin==3
        [dim,opts] = deal(opts,dim);
    else
        [dim,opts] = deal([],dim);
    end
end

% Row and column size {row,col} for the requested option.
switch lower(char(opts))
    case 'tot'
        sz = {sum(P.dim_out),sum(P.dim_in)};
    case 'block'
        sz = {size(P.C,1),size(P.C,2)};
    case 'split'
        sz = {P.dim_out,P.dim_in};
    otherwise
        error('cdopvar:size:opts','Size option must be one of ''tot'', ''block'' or ''split''.')
end

% 'split' returns a cell for a single output without dim, since dim_out
% and dim_in need not have the same length.
if ~isempty(dim)
    varargout = {[sz{dim}]};
elseif nargout>=2
    varargout = sz;
elseif strcmpi(opts,'split')
    varargout = {sz};
else
    varargout = {[sz{:}]};
end

end
