function K = params2kernel(P,Zin,Zd,dumnames)
% K = PARAMS2KERNEL(P,ZIN,ZD,DUMNAMES) returns the kernel defined by the
% sparse coefficients P, as a 'polynomial' or 'dpvar' object in the dummy
% variables DUMNAMES. Inverse of KERNELS2PARAMS, for a single kernel.
%
% INPUTS
% - P:          m*(nd+1) x n*nZ sparse array of coefficients, as an element
%               of the field params of an 'intop' object;
% - Zin:        1 x N cell, with element i the nmons(i) x di array of
%               monomial degrees along variable i, as the field Z.in of an
%               'intop' object, so that nZ = prod(nmons);
% - Zd:         nd x 1 cellstr, the names of the decision variables;
% - dumnames:   1 x sum(di) cellstr, the names of the dummy variables,
%               ordered variable by variable, as in
%               [Kop.varnames.in_dum{:}];
%
% OUTPUTS
% - K:          m x n object of type 'dpvar' if the kernel has decision
%               variables, or 'polynomial' otherwise. It only involves the
%               monomials, variables, and decision variables that appear in
%               it. If there are no decision variables, and the kernel is
%               constant with more than 1e6 elements, K is the sparse m x n
%               array of type 'double' instead;
%
% NOTES
% - Only the nonzeros of P, and the monomials they involve, are used, so
%   that the cost does not depend on the number of decision variables nd,
%   nor on the number of monomials nZ, apart from the size of P.
% - A zero kernel is declared as the monomial Z.in{1}(1,:)*...*Z.in{N}(1,:)
%   with coefficient 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - params2kernel
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

N = numel(Zin);
nd = numel(Zd);
[nmons,strides,nZ] = Zdims(Zin);
m = size(P,1)/(nd+1);
n = size(P,2)/nZ;

[r,c,v] = find(P);
r = r(:);   c = c(:);   v = v(:);   % row vectors if P has one row
j = ceil(c/nZ);
z = c - (j-1)*nZ;
blk = ceil(r/(nd+1));
d = r - (blk-1)*(nd+1);             % 1: constant term, l+1: decision variable l

% Degrees of the monomials that are used in all dummy variables, with the
% first variable changing slowest in the Kronecker product basis
zused = unique(z);
is_zero = isempty(zused);
if is_zero
    zused = 1;
end
nzu = numel(zused);
dg = zeros(nzu,0);
for i=1:N
    zi = mod(floor((zused-1)/strides(i)),nmons(i)) + 1;
    dg = [dg,Zin{i}(zi,:)];
end

if nd==0 && (is_zero || ~any(dg(:))) && m*n>1e6
    % Constant kernel with many elements: keep sparse
    K = sparse(r,j,v,m,n);
    return
end
[~,pz] = ismember(z,zused);
is_used_var = any(dg,1);
dg = dg(:,is_used_var);
varname_k = reshape(dumnames(is_used_var),[],1);
[dused,~,id] = unique(d(d>1));
if isempty(dused)
    % Plain polynomial: monomial rows, columns (j-1)*m+i
    K = polynomial(sparse(pz,(j-1)*m+blk,v,nzu,m*n),dg,varname_k,[m,n]);
else
    nu = numel(dused);
    dnew = ones(size(d));       dnew(d>1) = 1 + id;
    C = sparse((blk-1)*(nu+1)+dnew,(j-1)*nzu+pz,v,m*(nu+1),n*nzu);
    K = dpvar(C,dg,varname_k,Zd(dused-1),[m,n]);
end

end
