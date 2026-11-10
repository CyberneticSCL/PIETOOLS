function Vx = quad2lin_v2(Pmat,ZopL,ZxL,ZopR,ZxR,opts)
% VX = QUAD2LIN_V2(PMAT,ZOPL,ZXL,ZOPR,ZXR,OPTS) converts a distributed
% polynomial in quadratic form,
%   V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>_{L2}
% to the linear form, computing a functional operator Cop such that
%   V(x) = Cop*unique(ZL(x)oZR(x))
%
% INPUTS
% - Pmat:   m x n 'double' or 'dpvar' object representing the coefficients
%           parameterizing V(x) in the quadratic format, or a cell array
%           with Pmat{i,j} the block coupling monomials i of ZxL and j of
%           ZxR;
% - ZopL:   m x d1 'tensopmat' object representing the coefficient
%           operators acting on the distributed monomials ZxL;
% - ZxL:    d1 x 1 'polyopvar' object representing a basis of distributed
%           monomials;
% - ZopR:   (optional) n x d2 'tensopmat' object representing the
%           coefficient operators acting on the distributed monomials
%           ZxR. If not specified, ZopR=ZopL and ZxR=ZxL, for a square
%           Pmat;
% - ZxR:    (optional) d2 x 1 'polyopvar' object representing a basis of
%           distributed monomials;
% - opts:   (optional) char. Pass 'skip_combine' to skip the per-block
%           combine_terms folding of repeated-state monomials. Defaults
%           to '';
%
% OUTPUTS
% - Vx:     1 x 1 'polyopvar' object representing the function
%               V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>_{L2}
%           in the linear format;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - quad2lin_v2
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
% CR, 09/07/2026: Initial coding - Routes vector-valued TDP terms through
%                   the local term linearizer, preserving their Kronecker
%                   output form. Main diff to quad2lin is that
%                   quad2lin_term_v2 is called instead of quad2lin.
% DJ, 10/05/2026: Add optional input to skip the 'combine_terms' call;
% DJ, 10/15/2026: Updated to use new polyopvar functions;

if nargin<6 || isempty(opts)
    opts = '';                                                              % DJ, 10/05/2026
end
skip_combine = strcmpi(opts,'skip_combine');

% Assume a symmetric operator if only one half is specified.
is_symmetric = false;
if (nargin<4 || isempty(ZopR)) && size(Pmat,1)==size(Pmat,2)
    ZopR = ZopL;
    ZxR = ZxL;
    is_symmetric = true;
elseif nargin<4 || isempty(ZopR)
    error("Insufficient input arguments.")
else
    [ZxL,ZxR] = common_states(ZxL,ZxR);
end

% Extract the monomials, and their numbers
nZL = size(ZxL.Z,1);
nZR = size(ZxR.Z,1);
nfL = num_factors(ZxL);
nfR = num_factors(ZxR);

% Get the dimensions of the different operators
blkdimL = zeros(nZL,1);
for ii=1:nZL
    tmp = ZopL.ops{ii,ii};
    blkdimL(ii) = size(tmp,1);
end
blkdimL_cum = [0;cumsum(blkdimL)];
blkdimR = zeros(nZR,1);
for ii=1:nZR
    tmp = ZopR.ops{ii,ii};
    blkdimR(ii) = size(tmp,1);
end
blkdimR_cum = [0;cumsum(blkdimR)];

% Initialize the distributed polynomial in the linear format, as the empty
% polynomial (adding to the double 0 would add a constant term)
Vx = [];
for ii=1:nZL
    % Extract the monomial operators acting on the ith distributed
    % monomials
    Zop_ii = monomial_poly(ZxL,ZxL.Z(ii,:),ZopL.ops(ii,ii));                % DJ, 10/14/2026
    for jj=(ii-1)*is_symmetric+1:nZR
        % Extract the monomial operators acting on the jth distributed
        % monomial
        Zop_jj = monomial_poly(ZxR,ZxR.Z(jj,:),ZopR.ops(jj,jj));
        % Extract the block of decision variables acting on monomials Zi
        % and Zj
        if isa(Pmat,'cell')
            Pij = Pmat{ii,jj};
        else
            Pij = Pmat(blkdimL_cum(ii)+1:blkdimL_cum(ii+1), blkdimR_cum(jj)+1:blkdimR_cum(jj+1));
        end
        % Convert to the linear format
        Kij = quad2lin_term_v2(Pij,Zop_ii,Zop_jj);                       % CR, 09/07/2026
        % Account for symmetry
        if is_symmetric && ii~=jj
            Kij = 2*Kij;
        end
        % The monomial of this term: the factors of the monomial in ZL,
        % followed by those in ZR (the order of the inputs of the operator Kij)
        % Keep the zeros of Z: the factor of a constant (degree 0) monomial is a multiplier
        % from R, num_factors(ZxL)(ii) factors in all, see `polyopvar`
        zL = [full(ZxL.Z(ii,:)),zeros(1,nfL(ii))];
        zR = [full(ZxR.Z(jj,:)),zeros(1,nfR(jj))];
        Zrow = [zL(1:nfL(ii)),zR(1:nfR(jj))];
        Kpoly_ij = monomial_poly(ZxL,sparse(Zrow),{Kij});                   % DJ, 10/14/2026
        % Get rid of duplicate terms in the functional, and add to Vx
        if ~skip_combine                                                   % DJ, 10/05/2026
            Kpoly_ij = combine_terms(Kpoly_ij);
        end
        if isempty(Vx)
            Vx = Kpoly_ij;
        else
            Vx = Vx + Kpoly_ij;
        end
    end
end

end


function P = monomial_poly(Zx,Zrow,Cops)                                    % DJ, 10/14/2026
% P = MONOMIAL_POLY(ZX,ZROW,COPS) returns the 'polyopvar' in the states of ZX,
% with the single monomial ZROW, and the operators COPS (1 x 1 cell of the
% operator acting on this monomial) as coefficient. The output space is set
% to that of the operators, as in `piesos_quadvar`.
P = Zx;
P.Z = Zrow;
P.C = tensopmat(Cops);
P.dims.out = P.C.dims.out;
P.vars.out = P.C.vars.out;
P.dom.out = P.C.dom.out;
P.depmat.out = P.C.depmat.out;
P.chkval = true;
end