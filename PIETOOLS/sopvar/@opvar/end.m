function out = end(P,k,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out = end(P,k,n) is used to establish the size of opvar object P
% associated with index k out of n (in indexing, e.g. P(i:end,j:end))
%
% INPUT
% P:    opvar class object
% k,n:  natural numbers
%
% OUTPUT
% out:  the number of rows (k=1) or columns (k=2) of the full operator P,
%       or their product for linear indexing (n=1)
%
% NOTES:
% - This counts rows and columns of the full operator, as P(I,J) and
%   size(P) do, whereas size(P,'block') returns the number of blocks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - end (opvar)
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/14/2026: Initial coding;

% Rows times columns for a linear index, else rows (k=1) or columns (k=2)
if n==1
    out = sum(P.dim_out)*sum(P.dim_in);
elseif n==2 && k==1
    out = sum(P.dim_out);
elseif n==2
    out = sum(P.dim_in);
else
    error('opvar:endIndex','An opvar is indexed as P(I,J), or P(I) for one row or column.')
end

end
