function [Top_x,Fop,Kop,info] = compute_PIETOOLS_PDE_periodic_map_2D(PDE)
% Reconstruct a 2D PDE that has generalized periodic BCs along one of its
% two spatial directions, and regular (invertible) BCs along the other.
% We call p the periodic direction and q the other direction. Along p, the
% BCs are singular, and we introduce moment states v0 = Fop*x, which are
% functions of q only, following [Jagt, 2025]. We return operators such that
%   x = Top_x*[Dq*v0; xf],      Kop*[Dq*v0; xf] = 0,
% where xf is the highest-order mixed derivative of the PDE state, and Dq*v0
% is the highest-order derivative of v0 with respect to q.
%
% INPUTS:
% - PDE:    Initialized 'pde_struct' with two spatial variables, in which
%           all state components depend on both variables, and with BCs
%           not yet expanded (see "expand_PIETOOLS_PDE_BCs"), none of which
%           involves an input.
%
% OUTPUTS:
% - Top_x:  'opvar2d' mapping [Dq*v0; xf] to the PDE state x;
% - Fop:    'opvar2d' mapping the PDE state x to the moment states v0, which
%           are integrals of x along p;
% - Kop:    'opvar2d' defining the constraint Kop*[Dq*v0; xf]=0, which takes
%           values in L2[q];
% - info:   Struct with fields reconstruction (the operator Top_x), n_aux,
%           boundary_rank, jordan_transform, free_boundary_columns,
%           source_components, as in "compute_PIETOOLS_PDE_periodic_map",
%           and additionally
%               direction:      index (1 or 2) of the periodic variable;
%               selector:       'opvar2d' mapping [Dq*v0; xf] to v0 (or, if
%                               periodic in both directions, mapping the
%                               PIE state to all moment states);
%               x_tab_rows:     table rows describing the moment states,
%                               in the format of PDE.x_tab;
%
% NOTES:
% - The map is obtained by composing two 1D expansions. Along p, we expand
%       x = Tp*[v0; Dp*x],     K1*[v0; Dp*x] = 0,
%   using the 1D periodic converter, for each q. Along q, we expand
%       [v0; Dp*x] = Tq*[Dq*v0; Dq*Dp*x],
%   using the 1D converter with regular BCs. The latter requires the moment
%   states v0 to satisfy the same BCs as x. Since v0 = F*x with F acting
%   only along p, this is the case provided the coupling between the
%   different state components imposed by the BCs along q is also present
%   between the moment states, which we check.
% - Only BCs along p and q separately are supported (no mixed BCs), and
%   both inputs and ODE states are not supported.
% - If the BCs along both directions are singular, we obtain the moment
%   states v0 = F1*x (a function of s2), b = F2*x (a function of s1) and
%   z = F2*F1*x (finite-dimensional), and the PIE state is
%       [z; D1*b; D2*v0; xf].
%   This requires the BCs to be decoupled between the state components, see
%   "periodic_map_both_2D".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - compute_PIETOOLS_PDE_periodic_map_2D
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
% Initial coding: Claude 10/28/2026
% Claude, 10/29/2026: Add generalized PBCs along both spatial directions

id_unsup = 'PIETOOLS:convert:Periodic2D';
dom = PDE.dom;      vars = PDE.vars;
x_tab = PDE.x_tab;
if size(vars,1)~=2
    error(id_unsup,'Periodic boundary conditions are only supported in 2D by this routine.')
end
if ~all(x_tab(:,3) & x_tab(:,4))
    error(id_unsup,'Periodic boundary conditions in 2D are only supported if all state components depend on both spatial variables.')
end
ncomps = size(x_tab,1);     n = sum(x_tab(:,2));
for ii=1:numel(PDE.BC)
    for jj=1:numel(PDE.BC{ii}.term)
        if ~isfield(PDE.BC{ii}.term{jj},'x')
            error(id_unsup,'Periodic boundary conditions in 2D are not supported with inputs in the boundary conditions.')
        end
    end
end

%% Expand the state along each direction separately, using the 1D converter
PIE_dir = cell(1,2);    is_per = false(1,2);
for dd=1:2
    PDE_1D = build_1D_PDE(PDE,dd,id_unsup);
    PIE_dir{dd} = convert_PIETOOLS_PDE(PDE_1D,[],{'silent','Top'});
    is_per(dd) = ~isempty(PIE_dir{dd}.F);
end
if ~any(is_per)
    error('PIETOOLS:convert:NotPeriodic2D','The boundary conditions are not periodic along either spatial direction.')
elseif all(is_per)
    % Periodic along both directions
    [Top_x,Fop,Kop,info] = periodic_map_both_2D(PDE,PIE_dir,id_unsup);
    return
end
p = find(is_per);       q = 3-p;
PIEp = PIE_dir{p};      PIEq = PIE_dir{q};
Tp1 = PIEp.T;           Fp1 = PIEp.F;           Kp1 = PIEp.K;
info1 = PIEp.misc.periodic;
m = info1.n_aux;
vars_old_p = Tp1.var1;
vars_old_p = [Tp1.var1; Tp1.var2];
vars_old_q = [PIEq.T.var1; PIEq.T.var2];
vars_new_p = vars(p,:)';        vars_new_q = vars(q,:)';

%% Which scalar state does each moment state average, and its order in q
sizes = x_tab(:,2);
comp_of = repelem((1:ncomps)',sizes);
dq_scalar = x_tab(comp_of,4+q);
Fker = Fp1.Q1;
sigma = zeros(m,1);
for jj=1:m
    nz = find(~arrayfun(@(c) is_zero_poly(Fker(jj,c)),1:n));
    if numel(nz)~=1
        error(id_unsup,'Unexpected structure of the moment functional; each moment must average a single state component.')
    end
    sigma(jj) = nz;
end

%% Spaces and dimensions
% Block index of the space L2[q], in the ordering [R, L2[s1], L2[s2], L2[s1,s2]]
qb = 1+q;
dimMom = zeros(4,1);    dimMom(qb) = m;
dimX = [0;0;0;n];
kdim = Kp1.dim(1,1);
dimK = zeros(4,1);      dimK(qb) = kdim;
dimTop = dimMom + dimX;

%% Map along p:   x = Tp*[v0; w],  with w=Dp*x
Tp = opvar2d([],[dimX,dimTop],dom,vars);
if m>0
    Tp = set_param(Tp,4,qb,1,lift(Tp1.Q2,vars_old_p,vars_new_p));
end
Rp = Tp1.R;
for ii=1:3
    idx = [1,1];    idx(p) = ii;
    Tp.R22{idx(1),idx(2)} = lift(Rp.(['R',num2str(ii-1)]),vars_old_p,vars_new_p);
end

%% Map along q:   [v0; w] = Tq*[Dq*v0; Dq*w]
Rq = PIEq.T.R;
Rq_cell = cell(1,3);
for ii=1:3
    Rq_cell{ii} = polynomial(Rq.(['R',num2str(ii-1)]));
end
% The moment states satisfy the BCs along q in the same way as the
% corresponding components of x, if the BCs couple only moments that are
% averages with the same functional.
M_cell = cell(1,3);
for ii=1:3,     M_cell{ii} = polynomial(zeros(m,m));    end
for jj=1:m
    for kk=1:n
        coupled = false;
        for ii=1:3
            coupled = coupled || ~is_zero_poly(Rq_cell{ii}(sigma(jj),kk));
        end
        if ~coupled,    continue,   end
        if kk==sigma(jj)
            jk = jj;
        else
            % Find a moment that averages component kk using the same functional
            jk = 0;
            for ll=1:m
                if sigma(ll)==kk && is_zero_poly(Fker(ll,kk)-Fker(jj,sigma(jj)))
                    jk = ll;    break
                end
            end
            if jk==0
                error(id_unsup,['The boundary conditions along the non-periodic direction couple state components',...
                    ' that cannot be treated with the periodic boundary conditions along the other direction.'])
            end
        end
        for ii=1:3
            M_cell{ii}(jj,jk) = Rq_cell{ii}(sigma(jj),kk);
        end
    end
end
Tq = opvar2d([],[dimTop,dimTop],dom,vars);
Sel = opvar2d([],[dimMom,dimTop],dom,vars);
for ii=1:3
    Tq = set_param(Tq,qb,qb,ii,lift(M_cell{ii},vars_old_q,vars_new_q));
    Sel = set_param(Sel,qb,qb,ii,lift(M_cell{ii},vars_old_q,vars_new_q));
    idx = [1,1];    idx(q) = ii;
    Tq.R22{idx(1),idx(2)} = lift(Rq_cell{ii},vars_old_q,vars_new_q);
end

%% Compose
Top_x = Tp*Tq;

Fop = opvar2d([],[dimMom,dimX],dom,vars);
Fop = set_param(Fop,qb,4,1,lift(Fp1.Q1,vars_old_p,vars_new_p));

Kfun = opvar2d([],[dimK,dimTop],dom,vars);
Kfun = set_param(Kfun,qb,4,1,lift(Kp1.Q1,vars_old_p,vars_new_p));
Kop = Kfun*Tq;

%% Collect the information
dep = [q==1,q==2];
x_tab_rows = zeros(m,6);
for jj=1:m
    diff_vals = zeros(1,2);     diff_vals(q) = dq_scalar(sigma(jj));
    x_tab_rows(jj,:) = [info1.source_components(jj),1,dep,diff_vals];
end
info = info1;
info.reconstruction = Top_x;
info.direction = p;
info.selector = Sel;
info.x_tab_rows = x_tab_rows;
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function PDE_1D = build_1D_PDE(PDE,var_idx,id_unsup)
% Build a 1D PDE along just the spatial direction var_idx, satisfying the
% same BCs along this direction as the 2D PDE, as in the construction of
% the 1D maps in "Compute_Tmap_PDE_2D_Separable". All components must
% depend on both variables, and no inputs are included.
alt_idx = 3-var_idx;
pde_struct PDE_1D;
PDE_1D.x = PDE.x;
for ii=1:numel(PDE.x)
    PDE_1D.x{ii}.dom = PDE.x{ii}.dom(var_idx,:);
    PDE_1D.x{ii}.vars = PDE.x{ii}.vars(var_idx,:);
    PDE_1D.x{ii}.diff = PDE.x{ii}.diff(var_idx);
end
PDE_1D.x_tab = PDE.x_tab(:,[1,2,2+var_idx,4+var_idx]);
% Only the differentiability matters, not the actual dynamics
for ii=1:numel(PDE_1D.x)
    PDE_1D.x{ii} = rmfield(PDE_1D.x{ii},'term');
    PDE_1D.x{ii}.term{1}.x = ii;
    PDE_1D.x{ii}.term{1}.D = PDE_1D.x{ii}.diff;
end
% Retain only the BCs along var_idx
BC_idcs = ~PDE.BC_tab(:,2+var_idx);
PDE_1D.BC = PDE.BC(BC_idcs);
for ii=1:numel(PDE_1D.BC)
    PDE_1D.BC{ii}.vars = [];    PDE_1D.BC{ii}.dom = [];
    for jj=1:numel(PDE_1D.BC{ii}.term)
        if PDE_1D.BC{ii}.term{jj}.D(alt_idx)
            error(id_unsup,'The boundary conditions of the specified PDE are not supported.')
        end
        PDE_1D.BC{ii}.term{jj}.D = PDE_1D.BC{ii}.term{jj}.D(var_idx);
        PDE_1D.BC{ii}.term{jj}.loc = PDE_1D.BC{ii}.term{jj}.loc(var_idx);
        PDE_1D.BC{ii}.term{jj}.I = PDE_1D.BC{ii}.term{jj}.I(var_idx);
    end
end
PDE_1D = initialize(PDE_1D,true);
PDE_1D = reorder_comps(PDE_1D,'all',true);
end


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function tf = is_zero_poly(P)
% Check whether a 'polynomial' is (numerically) zero.
if isa(P,'double')
    tf = all(P(:)==0);
    return
end
tf = isempty(P.coefficient) || all(abs(P.coefficient(:))<1e-12);
end


function val = lift(P,vars_old,vars_new)
% Express a 1D polynomial kernel in terms of the variables of the 2D PDE.
val = subs(polynomial(P),vars_old,vars_new);
end


function P = set_param(P,rb,cb,ii,val)
% Set the parameter of the 'opvar2d' P mapping block cb to block rb, with
% index ii along the shared variable, to the value val. Here, blocks 1 to 4
% correspond to the spaces R, L2[s1], L2[s2], L2[s1,s2].
names = {'0','x','y','2'};
fname = ['R',names{rb},names{cb}];
if iscell(P.(fname))
    P.(fname){ii} = val;
else
    P.(fname) = val;
end
end



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [Top_x,Fop,Kop,info] = periodic_map_both_2D(PDE,PIE_dir,id_unsup)
% Map for generalized periodic BCs along both spatial directions.
% With the 1D periodic expansions along each direction d,
%   x = T0_d*F_d*x + Tx_d*D_d*x,        K_d*[F_d*x; D_d*x] = 0,
% we define the moment states
%   a = F1*x  (function of s2),   b = F2*x  (function of s1),   z = F2*a.
% Expanding x along s1, then a and D1*x along s2, we find
%   x = T0_1*a + Tx_1*(T0_2*u + Tx_2*xf),    a = Q0*z + Q1*y,
%   b = P0*z + P1*u,
% in terms of the PIE state v = [z; u; y; xf], with u = D1*b = F2*D1*x,
% y = D2*a, and xf = D1*D2*x. This state must satisfy the constraints
%   K1_f*(T0_2*u + Tx_2*xf) = 0,    K2_f*xf = 0,    K2_f*y = 0,
% where K1_f and K2_f are the constraints from the 1D expansions.
% The expansion of a and b is only valid if the BCs do not couple the
% different state components, which we check.
dom = PDE.dom;      vars = PDE.vars;
x_tab = PDE.x_tab;
ncomps = size(x_tab,1);     n = sum(x_tab(:,2));
sizes = x_tab(:,2);
comp_of = repelem((1:ncomps)',sizes);

%% Extract the data of the 1D expansions
D = cell(1,2);
for dd=1:2
    PIEd = PIE_dir{dd};
    vars_old = PIEd.vars(1,:)';     vars_new = vars(dd,:)';
    Td = PIEd.T;
    L.F = lift(PIEd.F.Q1,vars_old,vars_new);
    L.T0 = lift(Td.Q2,vars_old,vars_new);
    L.K = lift(PIEd.K.Q1,vars_old,vars_new);
    L.R = cell(1,3);
    for ii=1:3
        L.R{ii} = lift(Td.R.(['R',num2str(ii-1)]),vars_old,vars_new);
    end
    L.m = PIEd.misc.periodic.n_aux;
    L.k = PIEd.K.dim(1,1);
    L.info = PIEd.misc.periodic;
    % Source scalar state of each moment, with its derivative order
    L.sig = zeros(L.m,1);
    for jj=1:L.m
        nz = find(~arrayfun(@(c) is_zero_poly(L.F(jj,c)),1:n));
        if numel(nz)~=1
            error(id_unsup,'Unexpected structure of the moment functional; each moment must average a single state component.')
        end
        L.sig(jj) = nz;
    end
    % Check that the BCs do not couple different scalar states
    for ss=1:n
        for kk=1:L.m
            if L.sig(kk)~=ss && ~is_zero_poly(L.T0(ss,kk))
                error(id_unsup,['Periodic boundary conditions along both spatial directions are only supported',...
                        ' if the boundary conditions are not coupled between the state components.'])
            end
        end
        for ss2=1:n
            if ss2~=ss
                for ii=1:3
                    if ~is_zero_poly(L.R{ii}(ss,ss2))
                        error(id_unsup,['Periodic boundary conditions along both spatial directions are only supported',...
                                ' if the boundary conditions are not coupled between the state components.'])
                    end
                end
            end
        end
    end
    for rr=1:L.k
        if nnz(~arrayfun(@(c) is_zero_poly(L.K(rr,c)),1:n))~=1
            error(id_unsup,['Periodic boundary conditions along both spatial directions are only supported',...
                    ' if the boundary conditions are not coupled between the state components.'])
        end
    end
    D{dd} = L;
end
d_ord = zeros(n,2);
for dd=1:2,     d_ord(:,dd) = x_tab(comp_of,4+dd);  end
m1 = D{1}.m;    m2 = D{2}.m;
sig1 = D{1}.sig;    sig2 = D{2}.sig;

%% Pairs (j,k) of moments along s1 and s2 that average the same scalar state
pairs = zeros(0,2);
for jj=1:m1
    for kk=1:m2
        if sig1(jj)==sig2(kk),  pairs(end+1,:) = [jj,kk];   end %#ok<AGROW>
    end
end
npr = size(pairs,1);

%% Operators in the form of 'opvar2d' objects
% Dimensions of the spaces R, L2[s1], L2[s2], L2[s1,s2]
dimX = [0;0;0;n];
dimCols = [npr;m2;m1;n];          % PIE state v=[z;u;y;xf]
Z = @(rd,cd) opvar2d([],[rd,cd],dom,vars);
% Moment states a=F1*x (in L2[s2]), b=F2*x (in L2[s1]) and z
% Lifts of the 1D operators along each direction
for dd=1:2
    qb = 1+(3-dd);
    vo = D{dd}.m;
    dimMom = zeros(4,1);    dimMom(qb) = vo;
    dimK = zeros(4,1);      dimK(qb) = D{dd}.k;
    T0 = Z(dimX,dimMom);
    if vo>0,    T0 = set_param(T0,4,qb,1,D{dd}.T0);    end
    Tx = Z(dimX,dimX);
    for ii=1:3
        idx = [1,1];    idx(dd) = ii;
        Tx.R22{idx(1),idx(2)} = D{dd}.R{ii};
    end
    Fd = Z(dimMom,dimX);
    if vo>0,    Fd = set_param(Fd,qb,4,1,D{dd}.F);      end
    Kd = Z(dimK,dimX);
    if D{dd}.k>0,   Kd = set_param(Kd,qb,4,1,D{dd}.K);  end
    LL(dd).T0 = T0;     LL(dd).Tx = Tx;     LL(dd).F = Fd;      LL(dd).K = Kd; %#ok<AGROW>
end

% Selectors of the components of the PIE state
Sel_u = Z([0;m2;0;0],dimCols);
if m2>0,    Sel_u = set_param(Sel_u,2,2,1,polynomial(eye(m2)));     end
Sel_x = Z(dimX,dimCols);
Sel_x.R22{1,1} = polynomial(eye(n));
Sel_z = Z([npr;0;0;0],dimCols);
if npr>0,   Sel_z.R00 = eye(npr);   end

% a = Q0*z + Q1*y
Q0 = polynomial(zeros(m1,npr));
for pp=1:npr
    jj = pairs(pp,1);   kk = pairs(pp,2);
    Q0(jj,pp) = D{2}.T0(sig1(jj),kk);
end
Sel_a = Z([0;0;m1;0],dimCols);
if m1>0
    if npr>0,   Sel_a = set_param(Sel_a,3,1,1,Q0);  end
    for ii=1:3
        Q1 = polynomial(zeros(m1,m1));
        for jj=1:m1,    Q1(jj,jj) = D{2}.R{ii}(sig1(jj),sig1(jj));  end
        Sel_a = set_param(Sel_a,3,3,ii,Q1);
    end
end
% b = P0*z + P1*u
P0 = polynomial(zeros(m2,npr));
for pp=1:npr
    jj = pairs(pp,1);   kk = pairs(pp,2);
    P0(kk,pp) = D{1}.T0(sig2(kk),jj);
end
Sel_b = Z([0;m2;0;0],dimCols);
if m2>0
    if npr>0,   Sel_b = set_param(Sel_b,2,1,1,P0);  end
    for ii=1:3
        P1 = polynomial(zeros(m2,m2));
        for kk=1:m2,    P1(kk,kk) = D{1}.R{ii}(sig2(kk),sig2(kk));  end
        Sel_b = set_param(Sel_b,2,2,ii,P1);
    end
end

% w1 = D1*x = T0_2*u + Tx_2*xf
W1 = LL(2).Tx*Sel_x;
if m2>0,    W1 = W1 + LL(2).T0*Sel_u;   end
% x = T0_1*a + Tx_1*w1
Top_x = LL(1).Tx*W1;
if m1>0,    Top_x = Top_x + LL(1).T0*Sel_a;     end

%% Moment operator x -> [z; b; a]
Fz = Z([npr;0;0;0],dimX);
if npr>0
    Fker = polynomial(zeros(npr,n));
    for pp=1:npr
        jj = pairs(pp,1);   kk = pairs(pp,2);   ss = sig1(jj);
        Fker(pp,ss) = D{1}.F(jj,ss)*D{2}.F(kk,ss);
    end
    Fz = set_param(Fz,1,4,1,Fker);
end
Fop = [Fz; LL(2).F; LL(1).F];

%% Constraints: K1_f*w1 = 0, K2_f*xf = 0, K2_f*y = 0
Kc1 = LL(1).K*W1;
Kc2 = LL(2).K*Sel_x;
% Constraints on y: for each row of K2_f and each moment along s1 averaging
% the same scalar state
K3rows = zeros(0,2);
for rr=1:D{2}.k
    ss = find(~arrayfun(@(c) is_zero_poly(D{2}.K(rr,c)),1:n));
    for jj=find(sig1(:)'==ss)
        K3rows(end+1,:) = [rr,jj]; %#ok<AGROW>
    end
end
nr3 = size(K3rows,1);
K3 = Z([nr3;0;0;0],dimCols);
if nr3>0
    K3ker = polynomial(zeros(nr3,m1));
    for ll=1:nr3
        rr = K3rows(ll,1);  jj = K3rows(ll,2);
        K3ker(ll,jj) = D{2}.K(rr,sig1(jj));
    end
    K3 = set_param(K3,1,3,1,K3ker);
end
Kop = [K3; Kc2; Kc1];

%% Collect the information
x_tab_rows = zeros(npr+m2+m1,6);
for pp=1:npr
    x_tab_rows(pp,:) = [comp_of(sig1(pairs(pp,1))),1,0,0,0,0];
end
for kk=1:m2
    x_tab_rows(npr+kk,:) = [comp_of(sig2(kk)),1,1,0,d_ord(sig2(kk),1),0];
end
for jj=1:m1
    x_tab_rows(npr+m2+jj,:) = [comp_of(sig1(jj)),1,0,1,0,d_ord(sig1(jj),2)];
end
info = struct();
info.reconstruction = Top_x;
info.n_aux = npr+m2+m1;
info.boundary_rank = [D{1}.info.boundary_rank, D{2}.info.boundary_rank];
info.jordan_transform = {D{1}.info.jordan_transform, D{2}.info.jordan_transform};
info.free_boundary_columns = {D{1}.info.free_boundary_columns, D{2}.info.free_boundary_columns};
info.source_components = x_tab_rows(:,1);
info.direction = [1,2];
info.selector = [Sel_z; Sel_b; Sel_a];
info.x_tab_rows = x_tab_rows;
end
