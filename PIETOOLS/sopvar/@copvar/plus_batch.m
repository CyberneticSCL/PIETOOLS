function C = plus_batch(varargin)
% C = PLUS_BATCH(P1,P2,...,PN) is P1 + P2 + ... + PN for 'copvar'
% containers on one block grid. Fixed containers carry no decision variable
% lists, so the chain has nothing to reconcile and is used as is. A call
% with any 'cdopvar' operand dispatches to '@cdopvar/plus_batch'
% (InferiorClasses), which reconciles the lists once.
%
% INPUTS
% varargin: N >= 1 'copvar' containers on the same block grid, mapping
%           between the same spaces with the same component dimensions (as
%           for '+');
%
% OUTPUTS
% C:        'copvar' equal to P1 + P2 + ... + PN;
%
% NOTES
% Dispatched to when every operand is a 'copvar'; with any 'cdopvar'
% operand see '@cdopvar/plus_batch'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus_batch
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
% Initial coding MMP, 10/01/2026

if ~all(cellfun(@(x) isa(x,'copvar'),varargin))
    error('plus_batch:badInput','Summands must be copvar or cdopvar objects.')
end
C = varargin{1};
for k = 2:nargin,   C = C + varargin{k};    end
end
