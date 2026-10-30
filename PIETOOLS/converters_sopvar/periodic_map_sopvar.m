function [Top, Fop, Kop, info] = periodic_map_sopvar(PDE,Pop_f2x,Pop_b2x,bc_state_tab,Pop_b2BC,Pop_f2BC)
% Reconstruct a 1D PDE with singular homogeneous boundary conditions and
% mixed spatial orders, using 'copvar' operators.
% The FTC expansion gives x = Pop_b2x*xb + Pop_f2x*xf and BCs
%   0 = Pop_b2BC*xb + Pop_f2BC*xf,
% where xf takes the declared derivative of each state, and Pop_b2BC is a
% singular matrix. We retain additional moment states v0 = Fop*x, and return
% operators such that
%   x = Top.x*[v0;xf],      Kop*[v0;xf] = 0.
%
% INPUTS:
% - PDE:            Initialized 'pde_struct' with one spatial variable, and
%                   no ODE states;
% - Pop_f2x, Pop_b2x, bc_state_tab:  Expansion of the PDE state in terms of
%                   the fundamental state and core boundary state, see
%                   "ftc_expansion_sopvar";
% - Pop_b2BC, Pop_f2BC:     'copvar' operators describing the BCs, see
%                   "bc_operators_sopvar", with Pop_b2BC a singular matrix;
%
% OUTPUTS:
% - Top:    Struct with 'copvar' field x, mapping [v0;xf] to the PDE state,
%           and zero operators u and w;
% - Fop:    'copvar' mapping the PDE state to the moment states v0;
% - Kop:    'copvar' defining the constraint Kop*[v0;xf]=0;
% - info:   Struct with fields reconstruction (the operator Top.x), n_aux,
%           boundary_rank, jordan_transform, free_boundary_columns and
%           source_components.
%
% NOTES:
% - This implements Definition 1, Lemma 5, Proposition 6 and Lemma 7 of
%   Jagt, Chernyshenko and Peet, arXiv:2503.22896v5, as in
%   "compute_PIETOOLS_PDE_periodic_map", expressing all compositions of
%   operators using 'copvar' objects. The matrix H in the paper uses the
%   opposite sign to that of Pop_f2BC.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - periodic_map_sopvar
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
% Initial coding: Claude 10/28/2026, adapting compute_PIETOOLS_PDE_periodic_map to 'copvar'

vars = PDE.vars;        dom = PDE.dom;
s = vars(1,1);          theta = vars(1,2);
a = dom(1,1);           b = dom(1,2);           L = b-a;
x_tab = PDE.x_tab;
n = sum(x_tab(:,2));
G = full(copvar2matrix_local(Pop_b2BC));
nb = size(G,2);

%% Taylor monomials spanned by the core boundary state, Q(s)
% Each core boundary state is a derivative of a PDE state component at the
% lower boundary, with the PDE state recovered as (s-a)^k/k! times it.
nnx = [0;cumsum(x_tab(:,2))];
Q = polynomial(zeros(n,nb));
col_scalar = zeros(nb,1);       col_comp = zeros(nb,1);     col_order = zeros(nb,1);
cc = 0;
for kk=1:size(bc_state_tab,1)
    comp = bc_state_tab(kk,1);      nxc = bc_state_tab(kk,2);
    kder = bc_state_tab(kk,4);      % derivative order wrt the lower-boundary variable
    for t=1:nxc
        cc = cc+1;
        i = find(x_tab(:,1)==comp,1);
        col_scalar(cc) = nnx(i)+t;      col_comp(cc) = i;       col_order(cc) = kder;
        Q(col_scalar(cc),cc) = (s-a)^kder/factorial(kder);
    end
end

%% Lemma 3: Gauss-Jordan elimination with partial row pivoting
J = eye(size(G,1));
tol = eps;
pivots = zeros(1,0);    r = 0;
for j=1:nb
    [pivot,idx] = max(abs(G(r+1:end,j)));
    if isempty(pivot) || pivot<=tol
        continue
    end
    idx = idx+r;
    r = r+1;
    G([r,idx],:) = G([idx,r],:);
    J([r,idx],:) = J([idx,r],:);
    pivot = G(r,j);
    G(r,:) = G(r,:)/pivot;
    J(r,:) = J(r,:)/pivot;
    factors = G(:,j);   factors(r) = 0;
    G = G-factors*G(r,:);
    J = J-factors*J(r,:);
    pivots(end+1) = j; %#ok<AGROW>
    if r==nb,   break;  end
end
% Nonpivot columns identify the core-boundary states
free = setdiff(1:nb,pivots,'stable');
m = numel(free);
if m==0
    error('Boundary rank is numerically ambiguous. Check the boundary conditions.')
end
G1 = G(1:r,:);

%% Lemma 5: choose m auxiliary average states
source_scalars = col_scalar(free);
source_components = col_comp(free);
trace_orders = col_order(free);
is_second_order = all(x_tab(:,4)==2);
if is_second_order
    % Values of the Taylor monomials and their derivatives at the lower boundary
    Vm = zeros(n,nb);   Dm = zeros(n,nb);
    for j=1:nb
        Vm(col_scalar(j),j) = double(col_order(j)==0);
        Dm(col_scalar(j),j) = double(col_order(j)==1);
    end
    F3 = Vm(:,free)'/L + 2*(s-a)*Dm(:,free)'/L^2;
else
    F3 = Q(:,free)'/L;
end
G3 = double(int(F3*Q,s,a,b));
if rank([G1;G3])<nb
    if is_second_order
        F3 = 2*(2*b+a-3*s)*Vm(:,free)'/L^2 ...
            -6*(b+a-2*s)*Dm(:,free)'/L^3;
    else
        % Use polynomial moments dual to the free endpoint traces,
        % expanded in shifted Legendre polynomials
        max_order = max(x_tab(:,4));
        P = cell(max_order,1);
        P{1} = polynomial(1);
        z = (2*s-a-b)/L;
        if max_order>1
            P{2} = z;
        end
        for ell=1:max_order-2
            P{ell+2} = ((2*ell+1)*z*P{ell+1}-ell*P{ell})/(ell+1);
        end
        F3 = polynomial(zeros(m,n));
        for j=1:m
            k = trace_orders(j);
            d = x_tab(source_components(j),4);
            for ell=k:d-1
                weight = (2*ell+1)*(-1)^(ell+k)*prod(ell-k+1:ell+k) ...
                    /(factorial(k)*L^(k+1));
                F3(j,source_scalars(j)) = F3(j,source_scalars(j)) ...
                    + weight*P{ell+1};
            end
        end
    end
    G3 = double(int(F3*Q,s,a,b));
end
Ghat = [G1;G3];
if rank(Ghat)<nb
    error('Unable to construct independent average states for the singular boundary conditions.')
end
Ginv = Ghat\eye(nb);

%% Compose the operators using 'copvar' objects
% Spaces in which the operators exist
x_masks = 1;        x_dims = n;
% Moment functional Fop: L2^n[s] -> R^m, with kernel F3
empty_e = struct('r',{},'c',{},'rows',{},'cols',{},'sub',{},'C',{});
ent = struct('r',0,'c',1,'rows',1:m,'cols',1:n,'sub',1,'C',subs(F3,s,theta));
Fop = assemble_copvar_blocks(ent,0,m,x_masks,x_dims,vars,dom);
H3 = Fop*Pop_f2x;                                   % L2 -> R^m
% x = Q*Ginv*[0;v0] + (Pf2x - Q*Ginv*[H1;H3]) * xf, where H1 contains the
% first r BCs after elimination
Gm = matR(Ginv(:,r+1:end),vars,dom);
core = Gm*H3;
if r>0
    H1 = matR(J(1:r,:),vars,dom)*Pop_f2BC;          % L2 -> R^r
    Gb = matR(Ginv(:,1:r),vars,dom);
    core = Gb*H1 + core;
end
Cop = -(Pop_b2x*core);
T0 = Pop_b2x*Gm;                                    % R^m -> L2
Tx1 = Pop_f2x + Cop;                                % L2 -> L2
Tx = [T0, Tx1];

% Constraint K*[v0;xf] = 0, with K = [0, H(r+1:end)]
Kfun = matR(J(r+1:end,:),vars,dom)*Pop_f2BC;
Zmm = assemble_copvar_blocks(empty_e,0,m,0,m,vars,dom);
Kop = [Zmm, Kfun];

Top = struct();
Top.x = Tx;
Top.u = [];     Top.w = [];

info = struct('reconstruction',Tx,'n_aux',m, ...
    'boundary_rank',r,'jordan_transform',J,'free_boundary_columns',free, ...
    'source_components',source_components);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function P = matR(M,vars,dom)
% 'copvar' operator R^size(M,2) -> R^size(M,1) defined by a matrix.
ent = struct('r',0,'c',0,'rows',1:size(M,1),'cols',1:size(M,2),'sub',ones(1,size(vars,1)),'C',polynomial(M));
P = assemble_copvar_blocks(ent,0,size(M,1),0,size(M,2),vars,dom);
end

function G = copvar2matrix_local(P)
% Matrix representing a container of operators between finite-dimensional spaces.
rr = [0;cumsum(P.dim_out(:))];      cc = [0;cumsum(P.dim_in(:))];
G = zeros(rr(end),cc(end));
for i=1:numel(P.dim_out)
    for j=1:numel(P.dim_in)
        blk = P.C{i,j};
        if isempty(blk),    continue,   end
        G(rr(i)+1:rr(i+1),cc(j)+1:cc(j+1)) = full(blk.params{1});
    end
end
end
