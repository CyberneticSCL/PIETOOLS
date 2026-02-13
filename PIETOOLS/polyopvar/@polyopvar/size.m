function varargout = size(f,dim)
% VARARGOUT = SIZE(F,DIM) returns the size of the matrix-valued distributed
% polynomial function defined by 'polyopvar' object F along dimension DIM
% 
% INPUTS:
% - f:      'polyopvar' object;
% - dim:    scalar value 1 or 2 (optional input):
%           - if 1, number of rows of f is returned as output;
%           - if 2, number of columns of f is returned as output;
%           - if not specified, full dimension of f will be returned;
% 
% OUTPUTS:
% out: scalar or 1x2 array
%       - if one output is called and dim is specified, the size of the
%           function along the specified dimension is returned;
%       - if one output is called and no dim is specified, a 1x2 array
%           specifying [row,col] dimension of the polynomial is returned;
%       - if two outputs are called, the first output will be assigned the
%           number of rows, and the second the number of columns;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - size
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
% DJ, 02/13/2026: Initial coding

if nargin ==1 % if no dim is specified, return both row and column numbers 
    out = f.matdim;
    if nargout==0 || nargout==1
        varargout{1}=out;
    elseif nargout==2   % allow row and column numbers to be output separately
        varargout{1} = out(1);
        varargout{2} = out(2);
    end 
elseif nargin == 2
    if dim==1 || dim==2
        out = f.matdim(dim);
    else
        error('Dimension must be 1 or 2');
    end
    varargout{1} = out;  
end

end