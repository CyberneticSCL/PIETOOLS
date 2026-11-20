function [nmons,strides,nZ] = Zdims(Zin)
% [NMONS,STRIDES,NZ] = ZDIMS(ZIN) returns, for the Kronecker product basis
% of monomials defined by the cell ZIN, the number of monomials NMONS(i)
% along each direction, the stride STRIDES(i) of each direction in the
% Kronecker product (first direction slowest), and the total number of
% monomials NZ = prod(NMONS).
%
% INPUTS
% - Zin:     1 x N cell, with element i a nmons(i) x deps(i) array, each
%            row specifying the degrees of one monomial along direction i;
%
% OUTPUTS
% - nmons:   1 x N array of integers, with nmons(i) the number of monomials
%            along direction i;
% - strides: 1 x N array of integers, with strides(i) the number of
%            monomials along the directions i+1,...,N, so that the index of
%            the monomial in the Kronecker product basis increases by
%            strides(i) if the monomial along direction i is incremented;
% - nZ:      scalar integer specifying the total number of monomials,
%            nZ = prod(nmons);
%
% NOTES
% - If N=0, nZ=1.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Zdims
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
% DJ, 10/15/2026: Initial coding;

N = numel(Zin);
nmons = zeros(1,N);
for i=1:N
    nmons(i) = size(Zin{i},1);
end
strides = ones(1,N);
for i=N-1:-1:1
    strides(i) = strides(i+1)*nmons(i+1);
end
nZ = prod(nmons);

end
