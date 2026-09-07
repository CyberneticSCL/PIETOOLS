function [Aout,Bout] = apply_basis_map(T,Ain,Bin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Aout,Bout] = apply_basis_map(T,Ain,Bin) re-expresses one affine
% coefficient  vec(C) = A + B'*d  on the merged monomial bases, using the
% map T returned by 'sync_basis'.
%
% INPUT
% T:        sparse map from 'sync_basis', or [] when it is the identity;
% Ain,Bin:  one parameter cell of an sdopvar, i.e. vec(C) and its q x nC
%           block of decision variable coefficients;
%
% OUTPUT
% Aout,Bout: the same coefficient on the merged bases;
%
% NOTES
% An empty T means the operand is already on the merged bases, in which case
% both multiplies are returned untouched rather than performed. That case is
% common -- accumulating blocks into one operator, or adding a promoted
% 'sopvar' built on the same bases -- and T is a full nC x nC map, so the
% no-op is not cheap if it is actually carried out.
%
% T acts on the vec index, which is the COLUMN index of B; the rows of B are
% the decision variables and are never touched, so the cost here does not
% scale with their number.
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
% Initial coding MMP, 09/07/2026

if isempty(T)
    Aout = Ain;     Bout = Bin;
else
    Aout = T*Ain;   Bout = Bin*T.';
end

end
