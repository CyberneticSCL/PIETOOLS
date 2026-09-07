function varargout = size(P,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m,n] = size(P) returns the dimensions of the operator represented by the
% sdopvar object P, i.e. P maps L_2^n[vars.in] to L_2^m[vars.out].
%
% INPUT
% P:    sdopvar class object
% dim:  (optional) 1 to request the number of rows, 2 the number of columns
%
% OUTPUT
% varargout: [m,n] = P.dims, or P.dims(dim) if dim is specified
%
% NOTES:
% Without this method 'size' reports the size of the MATLAB object array,
% which is [1,1] for a single operator regardless of its true dimensions.
% 'end' is built on 'size', so the two must be present together.
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

if ~isa(P,'sdopvar')
    error('To check validity input must be sdopvar object');
end

matdim = P.dims;

if nargout<=1 && nargin==1
    varargout = {matdim};
    return
elseif nargin==1
    dim = 1:2;
end
varargout{1} = matdim(1);
varargout{2} = matdim(2);
varargout = varargout(dim);

end
