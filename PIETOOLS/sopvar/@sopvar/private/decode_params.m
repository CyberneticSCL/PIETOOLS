function [i,j,v,degmat,names] = decode_params(P,k)
% [I,J,V,E,NAMES] = DECODE_PARAMS(P,K) decodes the nonzero coefficients of
% the parameter P.PARAMS{K} of the 'sopvar' object P into the terms
%   V(t)*Z_L(s)*Z_R(s')
% of the kernel of the operator, for use in 'params2kernel' and in the
% display of the operator.
%
% INPUTS
% - P:      'sopvar' object, with M output variables and N input variables;
% - k:      integer, the linear index of the parameter in P.params;
%
% OUTPUTS
% - i,j:    nt x 1 arrays, the row and column of the matrix-valued kernel
%           that term t belongs to;
% - v:      nt x 1 array, the coefficients of the terms;
% - degmat: nt x (M+N) array, the degrees of the terms in the M output
%           variables P.vars.out and in the N dummy variables of the input
%           variables P.vars.in. In a direction in which the parameter is a
%           multiplier, the degree in the dummy variable is added to that
%           in the output variable, as the multiplier identifies the two,
%           so that E is zero for that dummy variable;
% - names:  1 x (M+N) cellstr, the names of the variables: those in
%           P.vars.out, followed by those in P.vars.in with _dum appended;
%
% NOTES
% - The rows of the parameter are ordered with the matrix row outermost and
%   the monomial innermost, and the columns with the matrix column
%   outermost and the monomial innermost. The monomial basis is the
%   Kronecker product of the bases P.ZL (rows) and P.ZR (columns), with the
%   first variable changing slowest.
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

[r,c,v] = find(P.params{k});
r = r(:);   c = c(:);   v = v(:);       % row vectors if the parameter has one row
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
names = [vout,strcat(vin,'_dum')];

end
