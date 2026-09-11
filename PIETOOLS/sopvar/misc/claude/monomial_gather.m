function p = monomial_gather(Z,ord)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p = monomial_gather(Z,ord) gives the permutation of a tensor monomial index
% induced by reordering the directions of Z.
%
% INPUT
% Z:    1 x N cell of exponent columns, in the OLD direction order;
% ord:  1 x N, ord(k) = index in Z of the direction that moves to new
%       position k, i.e. Znew = Z(ord);
%
% OUTPUT
% p:    prod(cellfun(@numel,Z)) x 1 gather index: new monomial j is old
%       monomial p(j). So a coefficient's monomial axis is reordered by
%       C_new = C_old(p,:) once p is lifted over the matrix dimension.
%
% NOTES
% The monomial vector is kron(Z{1},...,Z{N}) with the FIRST direction
% outermost, so reordering the cells is a tensor transposition of that index.
% The permutation is obtained by building the full degree table both ways and
% matching rows, which is the construction 'UnionBasisMonomials' uses, rather
% than by stride arithmetic: the table lives on the small (monomial) axis,
% and a wrong stride here would be silent.
%
% Column k of the new table describes new position k, i.e. old direction
% ord(k); reordering the old table's columns the same way makes the two
% tables describe the same monomials with matching columns, differing only in
% row order. Column reordering does not disturb row order, so the matched row
% index is the old monomial index.
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
% Initial coding MMP, 09/09/2026: split out of 'canonical_var_order' so that
%                  '@sopvar/mtimes' can use it to align two operands stored
%                  in different canonical variable orders.

if numel(Z)<=1 || isequal(ord(:).',1:numel(Z))
    p = (1:prod([cellfun(@numel,Z),1])).';
    return
end

Dold = degree_table(Z);
Dnew = degree_table(Z(ord));

[tf,p] = ismember(Dnew,Dold(:,ord),'rows');
if ~all(tf)
    error("monomial_gather: reordering lost a monomial.")
end
p = p(:);

end


%%
function D = degree_table(Z)
% One row per monomial of kron(Z{1},...,Z{N}), one column per direction,
% first direction outermost. Same construction as 'UnionBasisMonomials'.

D = zeros(1,0);
for k = 1:numel(Z)
    D = [kron(D,ones(size(Z{k}))), kron(ones(size(D,1),1),Z{k})];           %#ok<AGROW>
end

end
