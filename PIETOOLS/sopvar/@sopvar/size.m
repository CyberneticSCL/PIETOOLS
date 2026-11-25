function varargout = size(P,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARARGOUT = SIZE(P,DIM) returns the matrix dimensions of the 'sopvar'
% object P, the numbers of output and input components [q,p] of the PI
% operator L_2^p -> L_2^q it represents: SIZE(P) returns [q,p], [Q,P] =
% SIZE(P) returns them separately, and SIZE(P,DIM) returns the entries DIM.
%
% INPUT
% - P:      'sopvar' object;
% - dim:    (optional) 1 or 2, or a vector of them, the dimensions to
%           return.
%
% OUTPUT
% - varargout: with one output and no DIM the 1 x 2 vector [q,p]; with two
%           outputs q and p; with DIM the requested entries of [q,p].
%
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
% Initial coding AT -- 01/28/2026
% AT, 11/25/2026: With one output, SIZE(P,DIM) returns all requested
%                 entries as one vector (SIZE(P,[1 2]) returned only q);
%                 an error for DIM other than 1, 2 or [1,2].

% The input must be a 'sopvar' object
if ~isa(P,'sopvar')
    error('To check validity input must be sopvar object');
end

% The dimensions [q,p]; the sums below are over a single row
matdim = P.dims;

% One output and no DIM: the vector [q,p]
if nargout<=1 && nargin==1
varargout = {sum(matdim,1)};
return
elseif nargin==1
dim = 1:2;
elseif isempty(dim) || ~isnumeric(dim) || ~all(ismember(dim(:),[1,2]))
    error('sopvar:size:badDim','Dimension must be 1, 2 or [1,2]: a sopvar has only two dimensions.');
end
% Separate outputs, restricted to DIM
varargout{1} = sum(matdim(:,1));
varargout{2} = sum(matdim(:,2));
if nargout<=1
    % One output: the requested entries as one vector, SIZE(P,[1 2]) = [q,p]
    varargout = {[varargout{dim}]};
else
    varargout = varargout(dim);
end
end
