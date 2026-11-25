function idx = end(PDE,k,n)
% IDX = END(PDE,K,N) returns the last index in position K of N indices
% into the 'pde_struct' object PDE, as in x(end), PDE(end) or PDE(i,end).
%
% INPUT
% - PDE:    'pde_struct' object being indexed;
% - k:      integer, the position of 'end' among the indices;
% - n:      integer, the number of indices;
%
% OUTPUT
% - idx:    integer, the last index in position k: for a PDE variable x
%           of size m indexed as x(end), idx = m; otherwise, matching the
%           (equation, term) indexing PDE(i,j) of subsref, the number of
%           equations or of terms;
%
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/24/2026: Initial coding;

% Components x(idx) of a PDE variable.
[is_var,obj] = is_pde_var(PDE);
if is_var && n==1
    idx = PDE.(obj){1}.size;
    return
end
% Equations and terms: PDE(j) is term j if there is only one equation,
% and equation j otherwise.
[n_eqs,n_terms] = size(PDE);
if (n==1 && n_eqs==1) || k>=2
    idx = n_terms;
else
    idx = n_eqs;
end

end
