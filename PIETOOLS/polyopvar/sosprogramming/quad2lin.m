function Vx = quad2lin(Pmat,ZopL,ZxL,ZopR,ZxR)
% VX = QUAD2LIN(PMAT,ZOP,ZX) converts a distributed polynomial in quadratic
% form,
%   V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>_{L2}
% to the linear form, computing a functional operator Cop such that
%   V(x) = Cop*unique(ZL(x)oZR(x))
%
% INPUTS
% - Pmat:   m x n 'double' or 'dpvar' object representing the coefficients
%           parameterizing V(x) in the quadratic format;
% - ZopL:   m x d1 'tensopvar' object representing the coefficient
%           operators acting on the distributed monomials ZxL;
% - ZxL:    d1 x 1 'polyopovar' object representing a basis of distributed
%           monomials;
% - ZopR:   n x d2 'tensopvar' object representing the coefficient
%           operators acting on the distributed monomials ZxR. Defaults to
%           ZopL if not specified;
% - ZxR:    d2 x 1 'polyopovar' object representing a basis of distributed
%           monomials;
%
% OUTPUTS
% - Vx:     1 x 1 'polyopvar' object representing the function
%               V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>_{L2}
%           in the linear format;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - quad2lin
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
% DJ, 01/21/2026: Initial coding


% Assume a symmetric operator if only one half is specified.
is_symmetric = false;
if nargin<=3 && size(Pmat,1)==size(Pmat,2)
    ZopR = ZopL;
    ZxR = ZxL;
    is_symmetric = true;
elseif nargin<=3
    error("Insufficient input arguments.")
else
    [ZxL,ZxR] = common_vars(ZxL,ZxR);
end

% Extract the degrees of the distributed monomials
degmatL = ZxL.degmat;
nZL = size(degmatL,1);
degmatR = ZxR.degmat;
nZR = size(degmatR,1);

% Compute the list of monomials that appear in the distributed polynomial.
degmat_full = repmat(degmatL,[nZR,1]) + kron(degmatR,ones(nZR,1));
[M,degmat_lin] = uniquerows_integerTable(degmat_full);  % degmat_full = M*degmat_lin
old2new_idcs = M*(1:size(degmat_lin,1))';               % degmat_full = degmat_lin(old2new_idcs,:);
%nZ_lin = size(degmat_lin,1);

% Get the dimensions of the different operators
blkdimL = zeros(nZL,1);
for ii=1:nZL
    tmp = ZopL(ii,ii);
    blkdimL(ii) = tmp.matdim(1);
end
blkdimL_cum = [0;cumsum(blkdimL)];
blkdimR = zeros(nZR,1);
for ii=1:nZR
    tmp = ZopR(ii,ii);
    blkdimR(ii) = tmp.matdim(1);
end
blkdimR_cum = [0;cumsum(blkdimR)];

% Initialize an empty polynomial in the same variables as ZxL
tmp_poly = ZxL;
tmp_poly.C = tensopvar();
tmp_poly.degmat = zeros(0,numel(ZxL.varname));
Vx = tmp_poly;
for ii=1:nZL
    % Extract the monomial operators acting on the ith distributed
    % monomials
    Zop_ii = ZxL;
    Zop_ii.degmat = degmatL(ii,:);
    Zop_ii.C = ZopL(ii,ii);
    for jj=(ii-1)*is_symmetric+1:nZR
        % Extract the monomial operators acting on the jth distributed
        % monomial
        Zop_jj = ZxR;
        Zop_jj.degmat = degmatL(jj,:);
        Zop_jj.C = ZopR(jj,jj);
        %Zop_jj = ZopR(jj,jj);
        % Extract the block of decision variables acting on monomials Zi
        % and Zj
        Pij = Pmat(blkdimL_cum(ii)+1:blkdimL_cum(ii+1), blkdimR_cum(jj)+1:blkdimR_cum(jj+1));
        % Convert to the linear format
        Kij = quad2lin_term(Pij,Zop_ii,Zop_jj);
        % Account for symmetry
        if is_symmetric && ii~=jj
            Kij.params = 2*Kij.params;
        end
        % Match the obtained kernels with the correct distributed monomial
        lidx = old2new_idcs((ii-1)*nZL+jj);
        Kpoly_ij = tmp_poly;
        Kpoly_ij.C.ops{1} = Kij;
        Kpoly_ij.degmat = degmat_lin(lidx,:);
        % Get rid of duplicate terms in the functional, and add to Vx
        Kpoly_ij = combine_terms(Kpoly_ij);
        Vx = Vx + Kpoly_ij;
    end
end

end