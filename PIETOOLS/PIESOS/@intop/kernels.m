function K = kernels(Kop)
% K = KERNELS(KOP) returns the kernels defining the functional KOP, as a
% cell of the same size as KOP.params, with element K{i1,...,iN} the m x n
% kernel K_{i1,...,iN} associated with KOP.params{i1,...,iN}, see the
% documentation of the 'intop' class.
%
% INPUTS
% - Kop:    m x n 'intop' object;
%
% OUTPUTS
% - K:      cell of the same size as Kop.params, with each element an m x n
%           object of type 'dpvar' if the kernel has decision variables, or
%           'polynomial' otherwise, in terms of the dummy variables of the
%           functional, Kop.varnames.in_dum. Each kernel only involves
%           those monomials and those decision variables that appear in it,
%           so that the decision variables of different kernels may differ;
%
% NOTES
% - If Kop has no decision variables, a constant kernel with more than 1e6
%   elements is returned as a sparse array of type 'double' instead.
% - The sparse coefficients Kop.params can be used directly, instead of the
%   'polynomial' and 'dpvar' objects.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kernels
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
% DJ, 10/15/2026: Initial coding;

% Check that the functional is properly specified
if ~Kop.chkval && ~isvalid(Kop)
    error("The functional is not properly declared.")
end
N = numel(Kop.vars.in);
vn = Kop.varnames.in_dum;
dumnames = [cell(1,0),vn{:}];
[nmons,~,nZ] = Zdims(Kop.Z.in);
nd = numel(Kop.Zd);
[m,n] = size(Kop);

% Degree matrix of the full basis, first variable slowest
degmat_full = zeros(1,0);
for i=1:N
    degmat_full = [repelem(degmat_full,nmons(i),1), repmat(Kop.Z.in{i},size(degmat_full,1),1)];
end

K = cell(size(Kop.params));
for k=1:numel(K)
    P = Kop.params{k};
    [r,c,v] = find(P);
    r = r(:);   c = c(:);   v = v(:);   % row vectors if P has one row
    j = ceil(c/nZ);
    z = c - (j-1)*nZ;
    blk = ceil(r/(nd+1));               % matrix row
    d = r - (blk-1)*(nd+1);             % 1: constant term, l+1: decision variable l
    zused = unique(z);
    if nd==0 && (isempty(zused) || ~any(any(degmat_full(zused,:)))) && m*n>1e6
        % Constant kernel with many elements: keep sparse
        K{k} = sparse(r,j,v,m,n);
        continue
    end
    if isempty(zused)
        zused = 1;
    end
    nzu = numel(zused);
    pz = zeros(nZ,1);       pz(zused) = 1:nzu;
    dg = degmat_full(zused,:);
    is_used_var = any(dg,1);
    dg = dg(:,is_used_var);
    varname_k = reshape(dumnames(is_used_var),[],1);
    % Decision variables with nonzero coefficients in this kernel
    [dused,~,id] = unique(d(d>1));
    if isempty(dused)
        % Plain polynomial: monomial rows, columns (j-1)*m+i
        K{k} = polynomial(sparse(pz(z),(j-1)*m+blk,v,nzu,m*n),dg,varname_k,[m,n]);
    else
        nu = numel(dused);
        dnew = ones(size(d));       dnew(d>1) = 1 + id;
        C = sparse((blk-1)*(nu+1)+dnew,(j-1)*nzu+pz(z),v,m*(nu+1),n*nzu);
        K{k} = dpvar(C,dg,varname_k,Kop.Zd(dused-1),[m,n]);
    end
end

end
