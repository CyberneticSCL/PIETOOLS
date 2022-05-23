function varargout = size(P,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varargout = size(P,dim) gives the size of PI operator P: R^p x L2^q to
% R^m x L2^n along dimensions dim
% 
% INPUT
% P: opvar class object
% dim: a scalar or 1-by-2 non-zero integer vector
% 
% OUTPUT
% varargout: one value if dim is scalar, two values if no dim is specified
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PIETools - size
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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

if ~isa(P,'opvar')
    error('To check validity input must be opvar object');
end

matdim = P.dim;


if nargout<=1 && nargin==1
varargout = {sum(matdim,1)};
return
elseif nargin==1
dim = 1:4;
end

varargout{1} = sum(matdim(:,1));
varargout{3} = sum(matdim(1,:));
varargout{2} = sum(matdim(:,2));
varargout{4} = sum(matdim(2,:));


varargout = varargout(dim);
end
