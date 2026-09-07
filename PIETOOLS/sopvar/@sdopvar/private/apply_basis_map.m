function [Aout,Bout] = apply_basis_map(T,Ain,Bin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Aout,Bout] = apply_basis_map(T,Ain,Bin) re-expresses one affine
% coefficient  vec(C) = A + B'*d  on the merged monomial bases, using the
% map T returned by 'sync_basis'.
%
% INPUT
% T:        the scatter from 'sync_basis' -- a struct with fields 'idx'
%           (destination vec position of each source position) and 'nC'
%           (length of the destination) -- or [] when it is the identity;
% Ain,Bin:  one parameter cell of an sdopvar, i.e. vec(C) and its q x nC
%           block of decision variable coefficients;
%
% OUTPUT
% Aout,Bout: the same coefficient on the merged bases;
%
% NOTES
% An empty T means the operand is already on the merged bases, in which case
% both arguments are returned untouched rather than remapped. That case is
% common -- accumulating blocks into one operator, or adding a promoted
% 'sopvar' built on the same bases.
%
% T acts on the vec index, which is the COLUMN index of B; the rows of B are
% the decision variables and are only carried through, so the cost here does
% not scale with their number.
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
    % T.idx is injective -- distinct monomials take distinct positions in
    % the merged basis -- so this is a scatter with no accumulation, and one
    % indexed assignment per output expresses it.
    %
    % Assignment rather than a triplet build: measured against both the old
    % sparse multiply and a sparse() over find(Bin) at q = 2000, nC 729 ->
    % 13824, over densities 1e-4 to 1. Assignment was fastest at every
    % point (0.05x to 0.47x of the multiply), while the triplet build lost
    % to the multiply once B was dense (1.28x at nnz 921k) because that
    % constructor must sort every entry.
    Aout = sparse(T.nC,1);              Aout(T.idx)   = Ain;
    Bout = sparse(size(Bin,1),T.nC);    Bout(:,T.idx) = Bin;
end

end
