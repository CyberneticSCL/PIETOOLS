function Top = tmap_separable_sopvar(PDE)
% Compute 'copvar' operators Top.x, Top.u and Top.w such that the PDE state
% x can be expressed in terms of the fundamental state xf and the inputs as
%   x = Top.x * xf + Top.u * u + Top.w * w;
% for a PDE with an arbitrary number of spatial variables and "separable"
% boundary conditions.
%
% INPUTS:
% - PDE:    An initialized 'pde_struct' with ordered state components, as
%           passed to the conversion routine after "reorder_comps".
%
% OUTPUTS:
% - Top:    Struct with 'copvar' fields x, u, w. The field u (w) is [] if
%           there are no inputs u (w). The fundamental state is the
%           highest-order mixed derivative of the PDE state, as in
%           "ftc_expansion_sopvar".
%
% NOTES:
% - This is the n-D, 'copvar' version of "Compute_Tmap_PDE_2D_Separable".
%   For each variable s_k, we consider the 1D PDE along s_k with those
%   boundary conditions that are imposed on the boundary of s_k, and
%   compute the map T_k such that
%       x = T_k*D_k*x + Tw_k*w + Tu_k*u,
%   using "tmap_sopvar", where D_k takes the highest-order derivative of x
%   wrt s_k. Since T_k acts only along s_k, and is the identity along the
%   other variables, we can combine
%       x = T_1*T_2*...*T_n * D_1*D_2*...*D_n x + ...
% - Each boundary condition must be a condition on the boundary of
%   precisely one variable, as for example
%       x(a,s2,s3)=0,   x_{s1}(b,s2,s3)+x(b,s2,s3)=0,   x(s1,c,s3)=0;
%   but not, for example, x(a,c,s3)=0, or x(a,s2,s3)+x(s1,c,s3)=0.
% - At present, inputs may only enter the boundary conditions along a
%   single spatial direction, to avoid having to take derivatives of the
%   inputs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - tmap_separable_sopvar
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
% Initial coding, adapting Compute_Tmap_PDE_2D_Separable to 'copvar' and an
% arbitrary number of variables, 10/28/2026

vars = PDE.vars;
nvars = size(vars,1);
x_tab = PDE.x_tab;
BC_tab = PDE.BC_tab;
if nvars==1
    % Without other variables, the map along the variable is the full map,
    % also allowing for generalized periodic boundary conditions.
    Top = tmap_sopvar(PDE);
    return
end
pow2 = 2.^(0:nvars-1)';

% % % Establish for each BC the direction along which it is a boundary condition
% For each term involving a state component x, the BC must be imposed at the
% boundary of each variable that x depends on, but that the BC does not.
BC_dir = zeros(numel(PDE.BC),1);
for ii=1:numel(PDE.BC)
    in_BC = logical(BC_tab(ii,3:2+nvars));
    bndry_vars = false(1,nvars);
    for jj=1:numel(PDE.BC{ii}.term)
        term = PDE.BC{ii}.term{jj};
        if isfield(term,'x')
            in_x = logical(x_tab(term.x,3:2+nvars));
            bndry_vars = bndry_vars | (in_x & ~in_BC);
        end
    end
    if nnz(bndry_vars)>1
        error('The input PDE has mixed boundary conditions, which is not supported.')
    elseif nnz(bndry_vars)==1
        BC_dir(ii) = find(bndry_vars);
    end
    % BCs that involve no boundary values of a state cannot be enforced
    % along any direction.
end

% % % Expand along each direction
Tx = cell(nvars,1);     Tu = cell(nvars,1);     Tw = cell(nvars,1);
for k=1:nvars
    [Tx{k},Tu{k},Tw{k}] = tmap_1D_along(PDE,k,find(BC_dir==k));
end

% % % Combine the expansions
% Inputs can only be dealt with without differentiation if they contribute
% to the map along only one direction, which we expand first.
has_input = false(1,nvars);
for k=1:nvars
    has_input(k) = ~(isempty(Tu{k}) || is_zero_copvar(Tu{k})) || ~(isempty(Tw{k}) || is_zero_copvar(Tw{k}));
end
if nnz(has_input)>1
    error('The provided boundary conditions require derivatives of the inputs to be taken; this is not currently supported.')
end
first = find(has_input,1);
if isempty(first),  first = 1;      end
order = [first, setdiff(1:nvars,first)];
Top = struct();
Top.x = Tx{order(1)};
for k=order(2:end)
    Top.x = Top.x*Tx{k};
end
Top.u = Tu{first};
Top.w = Tw{first};
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
function [Tx,Tu,Tw] = tmap_1D_along(PDE,k,BC_idcs)
% Compute the map from the fundamental state to the PDE state along the
% direction k, as 'copvar' operators in all the variables of the PDE.
vars = PDE.vars;        dom = PDE.dom;
nvars = size(vars,1);
x_tab = PDE.x_tab;      u_tab = PDE.u_tab;      w_tab = PDE.w_tab;
pow2 = 2.^(0:nvars-1)';

% % % Declare the 1D PDE along the direction k
PDE_1D = pde_struct();
PDE_1D.x = PDE.x;
for ii=1:numel(PDE.x)
    if x_tab(ii,2+k)
        p = sum(x_tab(ii,3:2+k));       % position of variable k among those of the component
        PDE_1D.x{ii}.dom = PDE.x{ii}.dom(p,:);
        PDE_1D.x{ii}.vars = PDE.x{ii}.vars(p,:);
        PDE_1D.x{ii}.diff = PDE.x{ii}.diff(p);
    else
        PDE_1D.x{ii}.dom = [];
        PDE_1D.x{ii}.vars = [];
        PDE_1D.x{ii}.diff = [];
    end
    PDE_1D.x{ii} = rmfield(PDE_1D.x{ii},'term');
    PDE_1D.x{ii}.term{1}.x = ii;
    PDE_1D.x{ii}.term{1}.D = PDE_1D.x{ii}.diff;
end
PDE_1D.w = PDE.w;
for ii=1:numel(PDE.w)
    if w_tab(ii,2+k)
        p = sum(w_tab(ii,3:2+k));
        PDE_1D.w{ii}.dom = PDE.w{ii}.dom(p,:);
        PDE_1D.w{ii}.vars = PDE.w{ii}.vars(p,:);
    else
        PDE_1D.w{ii}.dom = [];
        PDE_1D.w{ii}.vars = [];
    end
end
PDE_1D.u = PDE.u;
for ii=1:numel(PDE.u)
    if u_tab(ii,2+k)
        p = sum(u_tab(ii,3:2+k));
        PDE_1D.u{ii}.dom = PDE.u{ii}.dom(p,:);
        PDE_1D.u{ii}.vars = PDE.u{ii}.vars(p,:);
    else
        PDE_1D.u{ii}.dom = [];
        PDE_1D.u{ii}.vars = [];
    end
end

% % % Keep only the BCs along direction k, treating them as 1D conditions
PDE_1D.BC = PDE.BC(BC_idcs);
for ii=1:numel(PDE_1D.BC)
    PDE_1D.BC{ii}.vars = [];        PDE_1D.BC{ii}.dom = [];
    for jj=1:numel(PDE_1D.BC{ii}.term)
        term = PDE_1D.BC{ii}.term{jj};
        if isfield(term,'x')
            tab = x_tab;        comp = term.x;
        elseif isfield(term,'w')
            tab = w_tab;        comp = term.w;
        else
            tab = u_tab;        comp = term.u;
        end
        if tab(comp,2+k)
            p = sum(tab(comp,3:2+k));
            if isfield(term,'x')
                if any(term.D([1:p-1,p+1:end]))
                    error('The boundary conditions of the specified PDE are not supported.')
                end
                term.D = term.D(p);
            end
            term.loc = term.loc(p);
            term.I = term.I(p);
        else
            if isfield(term,'x'),   term.D = [];    end
            term.loc = [];
            term.I = {};
        end
        PDE_1D.BC{ii}.term{jj} = term;
    end
end
% Drop BCs that involve only states that do not depend on variable k
is_ode_BC = true(numel(PDE_1D.BC),1);
for ii=1:numel(PDE_1D.BC)
    for jj=1:numel(PDE_1D.BC{ii}.term)
        term = PDE_1D.BC{ii}.term{jj};
        if isfield(term,'x') && x_tab(term.x,2+k)
            is_ode_BC(ii) = false;
        end
    end
end
PDE_1D.BC = PDE_1D.BC(~is_ode_BC);

PDE_1D = initialize(PDE_1D,true);
try
    Top_1D = tmap_sopvar(PDE_1D);
    if isfield(Top_1D,'F')
        error('Generalized periodic boundary conditions are only supported for PDEs in one spatial variable.')
    end
catch ME
    error('The expansion of the PDE state in terms of the PIE state fails; perhaps the PDE is ill-posed, or not currently supported.\n%s',ME.message)
end

% % % Lift the 1D operators to operators in all the variables of the PDE
x_dep = logical(x_tab(:,3:2+nvars));
[x_masks,x_dims,x_cols] = space_layout(x_tab(:,2),x_dep,nvars);
[~,~,x_cols1] = space_layout(x_tab(:,2),x_dep(:,k),1);
ent = lift_operator(Top_1D.x,k,vars,x_dep,x_cols1,x_cols,x_dep,x_cols1,x_cols);
Tx = assemble_copvar_blocks(ent,x_masks,x_dims,x_masks,x_dims,vars,dom);
Tu = [];    Tw = [];
for nm = {'u','w'}
    if isempty(PDE.(nm{1})),    continue,   end
    tab = PDE.([nm{1},'_tab']);
    i_dep = logical(tab(:,3:2+nvars));
    [i_masks,i_dims,i_cols] = space_layout(tab(:,2),i_dep,nvars);
    [~,~,i_cols1] = space_layout(tab(:,2),i_dep(:,k),1);
    if isempty(Top_1D.(nm{1}))
        ent = struct('r',{},'c',{},'rows',{},'cols',{},'sub',{},'C',{});
    else
        ent = lift_operator(Top_1D.(nm{1}),k,vars,x_dep,x_cols1,x_cols,i_dep,i_cols1,i_cols);
    end
    T = assemble_copvar_blocks(ent,x_masks,x_dims,i_masks,i_dims,vars,dom);
    if strcmp(nm{1},'u'),   Tu = T;     else,   Tw = T;     end
end
end


function entries = lift_operator(Top1D,k,vars,r_dep,r_cols1,r_cols,c_dep,c_cols1,c_cols)
% Express the entries of the 1D operator Top1D along variable k as entries
% of an operator in all variables, with identity along all variables other
% than k.
%
% INPUTS:
% - Top1D:      'copvar' acting along variable k;
% - k:          Index of the variable along which Top1D acts;
% - vars:       nvars x 2 array of primary and dummy variables of the PDE;
% - r_dep, c_dep:  ncomps x nvars logical arrays, indicating which variables
%               the output and input components depend on;
% - r_cols1, c_cols1:  For each component, the indices within the 1D space
%               (the R or the L2[s_k] space) of the 1D operator;
% - r_cols, c_cols:    For each component, the indices within the space it
%               lives in, in the operator in all variables;
nvars = size(vars,1);
pow2 = 2.^(0:nvars-1)';
Pop = copvar2opvar(Top1D);
vars_old = [Pop.var1;Pop.var2];
vars_new = [vars(k,1);vars(k,2)];
entries = struct('r',{},'c',{},'rows',{},'cols',{},'sub',{},'C',{});
for i=1:size(r_dep,1)
    in_i = r_dep(i,:);
    for j=1:size(c_dep,1)
        in_j = c_dep(j,:);
        if in_i(k) && in_j(k)
            ri = r_cols1{i};    cj = c_cols1{j};
            Ks = {Pop.R.R0(ri,cj), Pop.R.R1(ri,cj), Pop.R.R2(ri,cj)};
        elseif in_i(k)
            Ks = {Pop.Q2(r_cols1{i},c_cols1{j})};
        elseif in_j(k)
            Ks = {Pop.Q1(r_cols1{i},c_cols1{j})};
        else
            Ks = {Pop.P(r_cols1{i},c_cols1{j})};
        end
        for kk=1:numel(Ks)
            Kmat = polynomial(Ks{kk});
            if ~any(any(Kmat.coefficient)),     continue,   end
            % No variable other than k may be integrated over
            other = true(1,nvars);      other(k) = false;
            if any(in_j & ~in_i & other)
                error('The boundary conditions of the specified PDE are not supported.')
            end
            % The kernel uses the primary variable for multipliers and
            % output-only variables, and the dummy variable for variables
            % that are only integrated over in the input
            if in_i(k) && in_j(k)
                Kmat = subs(Kmat,vars_old,vars_new);
            elseif in_i(k)
                Kmat = subs(Kmat,vars_old(1),vars_new(1));
            elseif in_j(k)
                Kmat = subs(Kmat,vars_old(1),vars_new(2));
            end
            sub = ones(1,nvars);
            if in_i(k) && in_j(k),  sub(k) = kk;    end
            entries(end+1) = struct('r',in_i*pow2,'c',in_j*pow2, ... %#ok<AGROW>
                        'rows',r_cols{i},'cols',c_cols{j},'sub',sub,'C',Kmat);
        end
    end
end
end


function tf = is_zero_copvar(P)
% True if all blocks of the container are structurally or numerically zero.
tf = true;
for k=1:numel(P.C)
    blk = P.C{k};
    if isempty(blk),    continue,   end
    for m=1:numel(blk.params)
        if any(any(blk.params{m}~=0)),  tf = false;     return,     end
    end
end
end
