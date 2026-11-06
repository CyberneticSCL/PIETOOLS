function [Tops,Kops,Fop,info,Qops] = solve_BVP_1D(PDE,diffs)
% [Tops,Kops,Fop,info,Qops] = SOLVE_BVP_1D(PDE,diffs) takes a 'pde_struct'
% PDE defining a 1D boundary value problem (BVP), such as one of the
% problems returned by split_BVPs, and computes the operators
% parameterizing all functions that satisfy the boundary conditions (BCs),
% written as
%   sum_{j=1}^{d} E_{i,j}*[D_{j}u_{j}(a); D_{j}u_{j}(b)]
%       + sum_{j=0}^{d} int_{a}^{b} F_{i,j}(s)*u_{j}(s) ds = w_{i},
% where u_{j} stacks the states of order j, D_{j}u = [u; u_{s}; ...;
% (d/ds)^{j-1}u], and w stacks the right-hand sides of the BCs. For all
% u=(u_{0},...,u_{d}) satisfying the BCs, with v0=Fop*u and v1=(u_{0};
% u_{1}'; ...; (d/ds)^{d}u_{d}), the fundamental state being [v0;v1],
%     u = Tops.T0*v0 + Tops.T1*v1 + Tops.Tw*w,
%     Kops.x*v1 + Kops.w*w = 0,
% and, for the derivatives requested in diffs and 0<=i<=j<=d,
%     (d/ds)^{i} u_{j} = Tops.R0{j+1,i+1}*v0 + Tops.R1{j+1,i+1}*v1
%                                               + Tops.Rw{j+1,i+1}*w.
% Conversely, for all v0 and v1 with Kops.x*v1+Kops.w*w=0, the function
% u = Tops.T0*v0+Tops.T1*v1+Tops.Tw*w satisfies the BCs, with Fop*u=v0 and
% (d/ds)^{j}u_{j}=v1_{j}: the fundamental state is free of all BCs except
% Kops.x*v1+Kops.w*w=0.
%
% INPUT
% - PDE:    'pde_struct' (or 'struct') object defining a BVP in a single
%           spatial variable s in [a,b], see also "@pde_struct/initialize".
%           All states should depend on s, and have order PDE.x{k}.diff.
%           The dynamics are ignored. Every boundary condition 0 = sum_j
%           term_j should consist of terms of the states only, of the form
%           - C*(d/ds)^{D}u_{k}(a)  or  C*(d/ds)^{D}u_{k}(b),  with
%             constant C and D less than the order of u_{k};
%           - int_{a}^{b} C(s)*u_{k}(s) ds, over the full domain;
%           and the number of conditions (rows) should be at least the
%           total order of the states, sum_k n_k*diff_k, the surplus
%           constraining the fundamental state. The right-hand side w of
%           the conditions is not specified in PDE, but is the input of the
%           operators Tops.Tw, Kops.w and Tops.Rw. It stacks the
%           right-hand sides of the conditions PDE.BC{1},PDE.BC{2},...,
%           e.g. the forcing returned by split_BVPs;
% - diffs:  (optional) the orders i of the derivatives (d/ds)^{i}u_{j} for
%           which to compute Tops.R0, Tops.R1 and Tops.Rw (and Qops.Rp):
%           [] (the default) for none, 'all' for 0,...,d, or an array of
%           integers in 0,...,d;
%
% OUTPUT
% All operators are 'copvar' objects, each a 1 x 1 container holding a
% single 'sopvar' block, so that every operator maps between just one input
% and one output space (R^{p} or L2^{p}[a,b]). The states are regrouped by
% order, u = [u_{0};u_{1};...;u_{d}], with u_{j} stacking the states of
% order j in the order in which they appear in PDE.x, see info.x_order.
% Here m = n_tot-rank(G) is the number of auxiliary values of the
% fundamental state (see info.m), mK = n_BC-rank(G) the number of constraints
% on it, n_BC the number of BCs, n_{j} the dimension of u_{j}, |n| = sum_j
% n_{j}, and n_tot = sum_j j*n_{j} the total order of the states.
% - Tops:   struct with fields
%               T0: R^{m} --> L2^{|n|}, the map v0 --> u (empty if m=0),
%               T1: L2^{|n|} --> L2^{|n|}, the map v1 --> u,
%               Tw: R^{n_BC} --> L2^{|n|}, the map w --> u,
%               R0, R1, Rw: (d+1) x (d+1) cells, with elements
%                   R0{j+1,i+1}: R^{m} --> L2^{n_{j}} (empty if m=0),
%                   R1{j+1,i+1}: L2^{|n|} --> L2^{n_{j}},
%                   Rw{j+1,i+1}: R^{n_BC} --> L2^{n_{j}},
%               mapping v0, v1 and w to (d/ds)^{i} u_{j}, for 0<=i<=j<=d.
%               Only the elements for the orders i in diffs are computed,
%               the others (and those with n_{j}=0 or i>j) are empty;
% - Kops:   struct with fields
%               x: L2^{|n|} --> R^{mK}, v1 --> int_{a}^{b} K(s)*v1(s) ds,
%               w: R^{n_BC} --> R^{mK}, w --> -J2*w,
%               u: R^{0} --> R^{mK}, zero (there are no inputs u),
%           such that Kops.x*v1+Kops.w*w=0 for the fundamental state of u,
%           as in solve_PIETOOLS_BVP. The constraint does not depend on v0.
%           All have no rows if mK=0;
% - Fop:    L2^{|n|} --> R^{m}, u --> v0 = int_{a}^{b} F3(s)*u(s) ds. Empty
%           if m=0. The choice of F3 is not unique: if m=|n|, F3 is the
%           mean of the state, F3 = I/(b-a), provided this defines a valid
%           map, and otherwise F3 = P'*Phi, see compute_Tmap_1D;
% - info:   struct with fields
%           - n:        1 x d+1 array, with n(j+1) the total dimension of
%                       the states of order j;
%           - x_order:  1 x nx array, with PDE.x{x_order(k)} the k-th state
%                       in the order of u, so that u = [PDE.x{x_order(1)};
%                       PDE.x{x_order(2)};...];
%           - E,F:      the matrix and polynomial defining the BCs as
%                       above;
%           - lab_E,lab_F:  1 x 2*n_tot and 1 x |n| arrays, giving for each
%                       column of E and F the index of the state in PDE.x
%                       that it acts on;
%           - m:        number of auxiliary values of the fundamental
%                       state, the dimension of Fop's output. No auxiliary
%                       values are needed if m=0;
% - Qops:   (optional) struct splitting the map w --> u through the Taylor
%           coefficients c=Cw*w at s=a of the polynomial part of u that is
%           due to w, so that Tops.Tw = Tp*Cw and Tops.Rw{j+1,i+1} =
%           Rp{j+1,i+1}*Cw. The entries of c are the boundary values
%           Lambda_a u = [D_{1}u_{1}(a);...;D_{d}u_{d}(a)] of single
%           derivatives of single components, instead of the arbitrary
%           combinations in the rows of w. Fields
%               Cw: n_tot x n_BC array of type 'double', with c = Cw*w,
%               Tp: R^{n_tot} --> L2^{|n|}, c --> ebar(s-a)*c, with ebar
%                   the Taylor basis of the differentiable components,
%               Rp: (d+1) x (d+1) cell, with Rp{j+1,i+1}: R^{n_tot} -->
%                   L2^{n_{j}} the map c --> (d/ds)^{i} u_{j} for the
%                   polynomial ebar(s-a)*c, computed for the orders i in
%                   diffs, and empty elsewhere. Only computed if requested;
%
% NOTES
% - The BCs are written as E*Lambda_a u + int F u = w for compute_Tmap_1D,
%   which converts them to canonical form and computes the operators. This
%   routine only checks that the BCs are of the form above. Whether the
%   BCs of several 1D problems combine to a well-posed ND problem, e.g.
%   that they are order compatible, is not checked.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - solve_BVP_1D
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
% DJ, 11/02/2026: Initial coding;
% DJ, 11/06/2026: Allow more BCs than the total order of the states. The
%   surplus BCs give constraints on the fundamental state (Kops.x), so that
%   the BC data E,F and the input w have n_BC rows (was n_tot);

% % % Make sure the BVP is initialized, and is a valid 1D BVP
PDE = pde_struct(PDE);
if ~PDE.is_initialized
    PDE = initialize(PDE,true);
end
if size(PDE.vars,1)~=1
    error("The BVP should depend on exactly one spatial variable.")
end
nx = numel(PDE.x);
if nargin<2
    diffs = [];
end
if nx==0 || ~all(PDE.x_tab(:,3))
    error("The BVP should involve at least one state, and all states should depend on the spatial variable.")
end
b = PDE.dom(2);
nk = PDE.x_tab(:,2);    dk = PDE.x_tab(:,4);        % state dimensions and orders
d = max(dk);
if d==0
    error("None of the states is differentiable, so that there are no boundary values to be constrained.")
end


% % % Regroup the states by their order (compute_Tmap_1D requires
% % % u=[u_0;...;u_d]), keeping the order within each group
[~,x_order] = sort(dk(:)');
n_vec = zeros(1,d+1);       pos = zeros(nx,1);
for k=1:nx
    pos(k) = n_vec(dk(k)+1);                % offset of u_k in u_{dk(k)}
    n_vec(dk(k)+1) = n_vec(dk(k)+1) + nk(k);
end
n_tot = sum((1:d).*n_vec(2:end));           % number of boundary values
n_all = sum(n_vec);                         % dimension of the state u
nBC = numel(PDE.BC);
nrow = cellfun(@(c) c.size,PDE.BC);         % number of rows in each BC
row0 = [0;cumsum(nrow(:))];
% Surplus BCs (row0(end)>n_tot) become constraints on the fundamental
% state, Kops.x*v1+Kops.w*w=0
if row0(end)<n_tot                                                          % DJ, 11/06/2026
    error(['The number of boundary conditions (',num2str(row0(end)),') should be at least the total order of the states (',num2str(n_tot),').'])
end


% % % Write the BCs as E*Lambda_a u + int_{a}^{b} F(s)*u(s) ds = w, with
% % % Lambda_a u = [D_1u_1(a);D_2u_2(a);...;D_1u_1(b);...] as in
% % % compute_canonical_BCs, so that E has 2*n_tot columns
s = PDE.vars(1,1);      sd = PDE.vars(1,2);
idx_E0 = cumsum([0,2*(1:d).*n_vec(2:end)]);     % start of block j in E
idx_F0 = cumsum([0,n_vec]);                     % start of block j in F
% One row of E and F for each of the row0(end)>=n_tot conditions
E = zeros(row0(end),2*n_tot);       F = polynomial(zeros(row0(end),n_all)); % DJ, 11/06/2026

% Label the columns of E and F by the state they act on
lab_E = zeros(1,2*n_tot);            lab_F = zeros(1,n_all);
for k=1:nx
    j = dk(k);      idx = pos(k)+(1:nk(k));
    lab_F(idx_F0(j+1)+idx) = k;
    for dd=0:j-1
        lab_E(idx_E0(j)+dd*n_vec(j+1)+idx) = k;
        lab_E(idx_E0(j)+(j+dd)*n_vec(j+1)+idx) = k;
    end
end

for kk=1:nBC
    rws = row0(kk)+(1:nrow(kk));
    Rsel = zeros(row0(end),nrow(kk));   Rsel(rws,:) = eye(nrow(kk));
    for jj=1:numel(PDE.BC{kk}.term)
        t = PDE.BC{kk}.term{jj};
        if ~isfield(t,'x')
            error("The boundary conditions should only involve the states, not inputs.")
        end
        k = t.x;    j = dk(k);      cols = pos(k)+(1:nk(k));
        C = polynomial(t.C);
        if isempty(t.I{1})
            % Boundary value of derivative t.D<j, in block j of
            % E = [D_ju_j(a); D_ju_j(b)], each of size j*n_vec(j+1)
            isb = double(polynomial(t.loc(1)))==b;
            c = idx_E0(j) + isb*j*n_vec(j+1) + t.D*n_vec(j+1) + cols;
            E(rws,c) = E(rws,c) + double(C);
        else
            % Integral of the state over [a,b]
            if t.D~=0
                error("A boundary condition involves an integral of a derivative of a state; this is not supported.")
            end
            if ismember(sd.varname{1},C.varname)
                C = subs(C,sd,s);
            end
            Csel = zeros(nk(k),n_all);      Csel(:,idx_F0(j+1)+cols) = eye(nk(k));
            F = F + Rsel*C*Csel;
        end
    end
end


% % % Compute the operators
orders = diffs;
if ischar(diffs) && strcmp(diffs,'all')
    orders = 0:d;
end
if ~isempty(orders) && (~isnumeric(orders) || any(orders(:)<0) || any(orders(:)>d) || any(orders(:)~=round(orders(:))))
    error('solve_BVP_1D:diffs','The orders of the derivatives should be integers between 0 and the highest order of the states.')
end
orders = unique(orders(:)');
if nargout>=5
    [Tops,Rops,Kops,Fop,Qops] = compute_Tmap_1D(n_vec,E,F,PDE.dom,[s,sd],orders);
else
    [Tops,Rops,Kops,Fop] = compute_Tmap_1D(n_vec,E,F,PDE.dom,[s,sd],orders);
end
% The maps to the derivatives are returned with those to u
Tops.R0 = Rops.R0;    Tops.R1 = Rops.R1;    Tops.Rw = Rops.Rw;

% The constraints Kx*v1 + Kw*w + Ku*u = 0, with Kw = -Kw of compute_Tmap_1D
cols = struct('x',struct('dims',n_all,'depmat',true),'w',struct('dims',row0(end),'depmat',false),...
               'u',struct('dims',zeros(0,1),'depmat',false(0,1)));
Kw = Kops.Kw;
vars_1D = {s.varname{1}};
empty_op = @(P,sp) zeros_copvar_sop(struct('out',P.dim_out,'in',sp.dims),struct('out',{arrayfun(@(i) P.vars(P.space_out(i,:)),1:numel(P.dim_out),'UniformOutput',false)},'in',{arrayfun(@(i) vars_1D(sp.depmat(i,:)),1:numel(sp.dims),'UniformOutput',false)}),struct('vars',{vars_1D},'dom',PDE.dom));
K1 = Kops.K1;
Kops = struct('x',K1,'w',[],'u',[]);
if isempty(K1)
    K1 = struct('dim_out',zeros(0,1),'vars',{cell(1,0)},'space_out',false(0,0));
    Kops.x = empty_op(K1,cols.x);
end
Kops.w = empty_op(K1,cols.w);
Kops.u = empty_op(K1,cols.u);
if ~isempty(Kw)
    Kops.w = -Kw;
end

m = 0;                                      % Fop is empty if m=0
if ~isempty(Fop)
    m = Fop.dim_out(1);
end
info = struct('n',n_vec,'x_order',x_order,'E',E,'F',F,'m',m,'lab_E',lab_E,'lab_F',lab_F);

end
