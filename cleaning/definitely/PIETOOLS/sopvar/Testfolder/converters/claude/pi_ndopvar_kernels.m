function Kcells = pi_ndopvar_kernels(P,d,vars_dum)
% KCELLS = PI_NDOPVAR_KERNELS(P,D,VARS_DUM) evaluates the kernels of an
% 'ndopvar' object at the given values of its decision variables, returning
% them as matrix-valued 'polynomial' objects.
%
% The kernels are read straight off the class definition,
%
%   R_j(s,t;d) = (Im o Zd(s))^T (Ik o [1;d])^T C_j (In o Zd_t(t)),
%
% with k = m*prod(deg+1), Zd(s) the complete monomial basis of degree deg in
% every variable, and Zd_t(t) the complete basis in the INTEGRAL directions
% only, since a multiplier direction carries no dummy variable monomials.
% Nothing here goes through a conversion routine, so this can be used as an
% independent reference for the ndopvar/sdopvar converters.
%
% INPUTS
% - P:          'ndopvar' object;
% - d:          numel(P.dvarname) x 1 array of decision variable values, in
%               the order of P.dvarname;
% - vars_dum:   (optional) 1 x N 'cellstr' naming the dummy variables to use
%               in the output. Defaults to the dummy variables of P. Supply
%               this to line the kernels up with those of another object;
%
% OUTPUTS
% - Kcells:     cell array of the same size as P.C, holding each parameter as
%               a P.dim(1) x P.dim(2) 'polynomial' in the primary variables
%               of P and the requested dummy variables;
%
% MMP, 08/29/2026: Initial coding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pi_ndopvar_kernels
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

if ~isa(P,'ndopvar') && ~isa(P,'nopvar')
    error("Input must be an 'ndopvar' or 'nopvar' object.")
end

N = numel(P.deg);
m = P.dim(1);       n = P.dim(2);
if isa(P,'ndopvar')
    q = numel(P.dvarname);
else
    q = 0;
end
d = reshape(d,[],1);
if numel(d)~=q
    error("Number of decision variable values must match numel(P.dvarname).")
end
w = [1;d];

vars_s = P.vars(:,1).varname;
vars_s = reshape(cellstr(string(vars_s)),1,[]);
if nargin<3 || isempty(vars_dum)
    vars_dum = P.vars(:,2).varname;
end
vars_dum = reshape(cellstr(string(vars_dum)),1,[]);

nZ = prod(P.deg+1);
sz_C = [size(P.C),1];

Kcells = cell(size(P.C));
for j = 1:numel(P.C)
    Cj = P.C{j};
    if isempty(Cj) || nnz(Cj)==0
        Kcells{j} = polynomial(zeros(m,n));
        continue
    end

    % Which directions carry an integral, and hence a dummy monomial
    idcs = cell(1,N);
    [idcs{:}] = ind2sub(sz_C,j);
    idcs = cell2mat(idcs);
    is_int = logical(idcs-1);
    nZt = prod(P.deg(is_int)+1);

    if ~isequal(size(Cj),[m*nZ*(q+1),n*nZt])
        error("P.C{"+num2str(j)+"} is "+mat2str(size(Cj))+", expected "...
              +mat2str([m*nZ*(q+1),n*nZt])+".")
    end

    % Contract the q+1 rows of each slot against [1;d], without densifying
    E = full(kron(speye(m*nZ),sparse(w.'))*Cj);

    ZS = full_monoms(P.deg,vars_s);
    ZT = full_monoms(P.deg(is_int),vars_dum(is_int));

    K = polynomial(zeros(m,n));
    for p = 1:m
        for r = 1:n
            K(p,r) = ZS.'*(E((p-1)*nZ+(1:nZ),(r-1)*nZt+(1:nZt))*ZT);
        end
    end
    Kcells{j} = K;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = full_monoms(degs,vnames)
% The complete monomial basis kron_i (v_i^0,...,v_i^deg_i), first variable
% the slowest index.

degs = reshape(degs,1,[]);
degmat = zeros(1,0);
for i = 1:numel(degs)
    col = (0:degs(i))';
    degmat = [kron(degmat,ones(numel(col),1)), kron(ones(size(degmat,1),1),col)];  %#ok<AGROW>
end
T = size(degmat,1);
if isempty(vnames)
    Z = polynomial(1);
else
    Z = polynomial(speye(T),degmat,reshape(cellstr(string(vnames)),[],1),[T,1]);
end

end
