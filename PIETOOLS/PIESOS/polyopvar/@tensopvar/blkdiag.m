function Pop = blkdiag(varargin)
% POP = BLKDIAG(VARARGIN) takes the diagonal concatenation of tensor-PI 
% operators listed in VARARGIN, e.g.
%   POP = [varargin{1},0; 0, blkdiag(varargin{2:end})];
%
% INPUTS
% - varargin:   'tensopvar' objects representing tensor-PI operators;
%
% OUTPUTS
% - Pop:        'tensopmat' object representing the diagonal concatenation
%               of the input operators;
%
% NOTES
% - The concatenation is taken by the 'tensopmat' version of blkdiag, so
%   inputs of type 'tensopmat', 'intop' or 'sopvar' are also supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - blkdiag
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

% Just call the 'tensopmat' version
Pop = blkdiag(tensopmat(varargin{1}),varargin{2:end});

end