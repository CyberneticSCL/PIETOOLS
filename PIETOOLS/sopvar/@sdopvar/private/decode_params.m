function [i,j,l,v,degmat] = decode_params(P,k)
% [I,J,L,V,DEGMAT] = DECODE_PARAMS(P,K) decodes the nonzero coefficients of
% the parameter K of the 'sdopvar' object P into the terms
%   V(t)*Zd(L(t))*Z_L(s)*Z_R(s')
% of the kernel of the operator, for use in the display of the operator.
%
% INPUTS
% - P:      'sdopvar' object, with M output variables and N input
%           variables;
% - k:      integer, the linear index of the parameter in P.params.A, with
%           coefficients vec(C) = P.params.A{k} + P.params.B{k}'*P.Zd;
%
% OUTPUTS
% - i,j:    nt x 1 arrays, the row and column of the matrix-valued kernel
%           that term t belongs to;
% - l:      nt x 1 array, the index in P.Zd of the decision variable that
%           term t is a multiple of, or 0 if term t is a constant;
% - v:      nt x 1 array, the coefficients of the terms;
% - degmat: nt x (M+N) array, the degrees of the terms in the M output
%           variables P.vars.out and in the N dummy variables of the input
%           variables P.vars.in. In a direction in which the parameter is a
%           multiplier, the degree in the dummy variable is added to that
%           in the output variable, as the multiplier identifies the two,
%           so that degmat is zero for that dummy variable;
%
% NOTES
% - The coefficient matrix C of the parameter has the matrix row outermost
%   and the monomial innermost in its rows, and the matrix column outermost
%   and the monomial innermost in its columns, with the monomials the
%   Kronecker product of the bases P.ZL (rows) and P.ZR (columns), the
%   first variable changing slowest. Only the nonzeros of P.params.A{k} and
%   P.params.B{k} are accessed, so the cost does not depend on the number
%   of decision variables.
% - The parameters are indexed by multi-indices over the variables shared
%   by the input and output, see 'gamma_of_cell': 1 for a multiplier, 2 for
%   the integral from the lower boundary to the variable, and 3 for the
%   integral from the variable to the upper boundary.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - decode_params
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
% DJ, 11/06/2026: Initial coding;

vout = reshape(P.vars.out,1,[]);
vin = reshape(P.vars.in,1,[]);
M = numel(vout);
N = numel(vin);
nL = cellfun(@numel,P.ZL);
nR = cellfun(@numel,P.ZR);
nZL = prod(nL);
nZR = prod(nR);
nrow = P.dims(1)*nZL;

% The constant terms, and those multiplying a decision variable, as indices
% into the vectorized coefficient matrix. A zero block may be stored as [] or 0
[idx_a,~,v_a] = find(P.params.A{k}(:));
[l_b,idx_b,v_b] = find(P.params.B{k});
idx = [idx_a(:);idx_b(:)];
l = [zeros(numel(idx_a),1);l_b(:)];
v = [v_a(:);v_b(:)];
r = rem(idx-1,nrow)+1;
c = ceil(idx/nrow);
i = ceil(r/nZL);
zl = r - (i-1)*nZL;
j = ceil(c/nZR);
zr = c - (j-1)*nZR;
aL = kron_split(zl-1,kron_strides(nL));
aR = kron_split(zr-1,kron_strides(nR));
degmat = zeros(numel(v),M+N);
for q=1:M
    degmat(:,q) = P.ZL{q}(aL(:,q)+1);
end
for q=1:N
    degmat(:,M+q) = P.ZR{q}(aR(:,q)+1);
end

% In a multiplier direction, add the degree in the dummy variable to that in
% the output variable
S3 = intersect(vin,vout);
if ~isempty(S3)
    gam = gamma_of_cell(k,numel(S3));
    for t=find(gam==1)
        pl = find(strcmp(vout,S3{t}));
        pr = find(strcmp(vin,S3{t}));
        degmat(:,pl) = degmat(:,pl) + degmat(:,M+pr);
        degmat(:,M+pr) = 0;
    end
end

end
