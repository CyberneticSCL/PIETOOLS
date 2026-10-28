function A = expand_full(a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A = EXPAND_FULL(a) returns the 2^nf multi-indices obtained from the
% multi-index a by replacing each of its nf entries 4 (a full-domain
% integral) by 2 (a lower integral) or by 3 (an upper integral), one per
% row, the first entry 4 varying fastest; A = a when a has no entry 4.
%
% INPUT
% - a:      'double', 1 x n, entries in {1,2,3,4}: per direction a
%           multiplier (1), a lower (2), an upper (3) or a full-domain (4)
%           integral.
%
% OUTPUT
% - A:      'double', 2^nf x n, entries in {1,2,3}.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - expand_full
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
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body; code unchanged. The deleted header text was
%                the reason a full-domain index stands for its 2/3
%                substitutions and a comparison with the layout of
%                'alpha_all'. The file had no log: initial coding MMP,
%                09/12/2026, as a local function of 'possopvar' (the file
%                since renamed 'sopquadvar'), moved to private/ by MMP,
%                09/21/2026.

% Directions with a full-domain integral; none: a itself
f = find(a==4);
if isempty(f)
    A = a;
    return
end
% One copy of a per substitution; column f(t) alternates 2 and 3 in runs
% of 2^(t-1) rows
nf = numel(f);
A = repmat(a,2^nf,1);
for t = 1:nf
    blk = 2^(t-1);
    A(:,f(t)) = repmat(reshape(repmat([2,3],blk,1),[],1),2^nf/(2*blk),1);
end

end
