function dV = Liediff(V,PIE)
% DV = LIEDIFF(V,PIE) returns the Lie derivative of a distributed
% polynomial functional V along the PIE defined by PIE,
%   d/dt (T*x(t)) = f(x(t)),
% so that dV(x) = d/dt V(x(t)) for solutions x of the PIE.
%
% INPUTS
% - V:      1 x 1 'polyopvar' object representing a distributed polynomial
%           functional (as declared by e.g. `piesos_quadvar`, `innerprod`),
%           in terms of the PDE state variables u=T*x, and potentially also
%           the associated PIE state variables x (see NOTES);
% - PIE:    'struct' representing a distributed polynomial PIE, as output
%           by `convert_PIETOOLS_PDE_nonlinear`, with fields
%           T:  'tensopmat' representing the map x -> u=T*x from PIE state
%               to PDE state;
%           f:  'polyopvar' representing the right-hand side of the PIE,
%               d/dt u = f(x);
%           x:  'polyopvar' array of the individual PIE state variables;
%           varname_pde, varname_pie: cellstr with the names of the PDE and
%               PIE state variables, in the same order as the rows of T, f;
%
% OUTPUTS
% - dV:     1 x 1 'polyopvar' object representing the derivative of V along
%           the PIE, in terms of the PIE state variables PIE.varname_pie
%           (V and dV as declared by `piesos_LF`);
%
% NOTES
% - Each monomial of V is treated separately:
%   - In PDE states only, V(u)=K(u_{i1} o ... o u_{id}): by the product
%     rule, dV(x) = sum_j K(u_{i1} o ... o f_{ij}(x) o ... o u_{id}), with
%     the remaining u replaced by T*x;
%   - In PDE and PIE states, V(x,u)=K(x_{i1} o ... o u_{j1} o ...): the
%     symmetry condition V(xL,T*xR)=V(xR,T*xL) for independent copies xL,
%     xR of the PIE state is assumed to hold, as enforced by `piesos_LF`
%     for PIE states. Then the derivative of the factors x, which are left
%     as is, equals that of the factors u, so that dV is twice the sum
%     over the factors u of the terms above. The symmetry is not verified:
%     without it, dV is not the derivative of V;
%   - In PIE states only: an error, as d/dt x is not defined by the PIE,
%     only d/dt (T*x). A constant term has zero derivative, and is
%     dropped.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Liediff
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

% Check the inputs
if ~isa(V,'polyopvar')
    error("The functional must be specified as 'polyopvar' object.")
end
if ~isequal(size(V),[1,1])
    error("The functional must be scalar-valued.")
end
if ~isstruct(PIE) || ~all(isfield(PIE,{'T','f','x','varname_pde','varname_pie'}))
    error("The PIE must be specified as 'struct' with fields 'T', 'f', 'x', 'varname_pde', and 'varname_pie'.")
end
if ~V.chkval
    V = initialize(V);          % guarantee V.Z and V.C are consistent
end
V = clean_poly(V);

% Every state in V must be a PDE state u or a PIE state x of the PIE; names
% in both lists could not be told apart
names_pde = PIE.varname_pde(:);
names_pie = PIE.varname_pie(:);
if any(ismember(names_pde,names_pie))
    error("PDE and PIE state variables must have different names.")
end
if ~all(ismember(V.varname,[names_pde;names_pie]))
    error("All state variables of the functional must appear as PDE or PIE state in the PIE.")
end

% Express the PDE states in V in terms of the PIE state, u=T*x, and declare
% their time derivatives, d/dt u=f(x)
idcs = find(ismember(names_pde,V.varname))';
if isempty(idcs)
    error("The functional does not depend on any PDE state. The derivative of a functional in PIE states only is not defined by the PIE, which specifies only d/dt (T*x).")
end
names_u = names_pde(idcs);
Top = PIE.T;
fx = PIE.f;
if ~isequal(idcs,1:numel(names_pde))
    % Restrict T, f to the rows of the PDE states in V, as in piesos_LF
    Top = Top(idcs,:);
    fx = fx(idcs,:);
end
Tx = Top*PIE.x;                     % row i is the PDE state names_u{i}
fu = cell(numel(idcs),1);           % fu{i}: d/dt of the PDE state names_u{i}
for i=1:numel(idcs)
    if size(fx,1)==1
        fu{i} = fx;                 % single row: use as is, as in piesos_LF
    else
        fu{i} = fx(i,:);
    end
end

% Take the derivative of each monomial of V
dVu = 0;                            % monomials in PDE states only
dVxu = 0;                           % monomials in PDE and PIE states
for k=1:size(V.Z,1)
    if isempty(V.C.ops{1,k})
        continue
    end
    Wk = V{1,k};
    zk = full(Wk.Z(1,:));
    pos = find(zk>0);               % positions of the state factors in Wk.Z (rest is padding)
    if isempty(pos)
        continue                    % constant: zero derivative
    end
    [is_u,rows] = ismember(Wk.varname(zk(pos)),names_u);
    if ~any(is_u)
        error("The derivative of a monomial in PIE states only is not defined by the PIE, which specifies only d/dt (T*x).")
    end
    % Product rule: replace each factor u in turn by f(x), and the other
    % factors u by T*x. Factors are numbered by position in Wk.Z, as in subs
    dWk = 0;
    for j=find(is_u(:))'
        dWk = dWk + subs(subs(Wk,pos(j),fu{rows(j)}),names_u,Tx);
    end
    if all(is_u)
        dVu = dVu + dWk;
    else
        dVxu = dVxu + dWk;          % factors x are not differentiated
    end
end
if ~isnumeric(dVxu)
    dVxu = 2*dVxu;                  % symmetry: the factors x give the same contribution
end
dV = dVu + dVxu;
if isnumeric(dV)
    error("The functional has no monomials in the PDE states, so that its derivative along the PIE is zero.")
end
dV = clean_poly(dV);

end
