function out = end(a,k,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = end(a,k,n) is used to establish the size of sdopvar object a
% associated with index k out of n (in indexing, e.g. a(i:end,j:end))
%
% INPUTS
% a:    'sdopvar' object
% k,n:  natural numbers, the position k of the index among n indices
%
% OUTPUTS
% out:  the size of a along dimension k (n = 2), or the number of entries,
%       rows times columns, for a linear index (n = 1)
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

% Number of entries for a linear index, else the size along dimension k
if n ==1
    sza = size(a);
    out = sza(1)*sza(2);
elseif n==2
    out = size(a,k);
else
    error('sdopvar:endIndex','An sdopvar is indexed as P(I,J), or P(I) for one row or column.')
end

end
