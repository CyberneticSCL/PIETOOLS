function dirs = split_BVPs(PDE)
% DIRS = SPLIT_BVPS(PDE) splits the boundary conditions (BCs) of a PDE
% into 1D problems, one along each spatial direction s_i, written as
%   E*[D_{j}u_{j}(a); D_{j}u_{j}(b)]_{j} + int_{a}^{b} F(s)*u(s) ds = g
% for compute_Tmap_1D. E.g. the BCs u(0,y)=u(1,y)=0, u(x,0)=u_{y}(x,1)=0
% of u(x,y) give u(0)=u(1)=0 along x, and u(0)=u_{y}(1)=0 along y. The
% dynamics are ignored.
%
% INPUT
% - PDE:    initialized 'pde_struct' with N spatial variables. Each BC
%           must be a condition along a single direction, see NOTES;
%
% OUTPUT
% - dirs:   1 x N struct array, describing the BCs along s_i on a slice
%           (the other variables fixed), with fields
%           - BC_idcs:  the BCs along s_i, with nrow their numbers of
%                       rows;
%           - x_idcs:   the states that depend on s_i;
%           - n:        n(j+1) the total size of the states of order j
%                       in s_i, with u = [u_0;...;u_d] stacking the
%                       states x_idcs(x_order);
%           - E, F:     the BCs as above, with lab_E, lab_F the state
%                       (index in x_idcs) of each column;
%           - forcing:  cell, with forcing{k} the terms of BC k that do
%                       not involve a state depending on s_i (e.g.
%                       inputs), with C negated: they sum to g;
%           If no state is differentiable in s_i (and there are no BCs
%           along s_i), n, x_order, E, F, lab_E and lab_F are empty;
%
% NOTES
% - A BC is a condition along s_i if it fixes s_i, by evaluating states
%   at a boundary or integrating over s_i, and varies in the other
%   variables of its states, e.g.
%       x(a,y)=0,       x_{x}(b,y)+int_{c}^{d} x(b,z)dz = w(y).
%   Not supported (error): BCs that fix several variables (corners or
%   edges; convert_PIETOOLS_BVP imposes these separately if each term of
%   a state does so) or that also vary in s_i, derivatives along the
%   boundary, coefficients of the states varying along the boundary,
%   states depending on s_i on different variables, and delays of
%   states. There should be at least as many conditions along s_i as the
%   total order of the states in s_i. Integrals of derivatives along s_i
%   are integrated by parts.
% - Whether the 1D problems combine to an admissible problem is checked
%   in compose_Tmaps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - split_BVPs
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
% DJ (CC), 11/02/2026: Initial coding;
% DJ (CC), 11/13/2026: Return E, F of the 1D problems (no 1D 'pde_struct');
%                   integrate integrals of derivatives by parts; shorter
%                   error message;

N = size(PDE.vars,1);
if N==0
    error("The PDE does not depend on any spatial variables.")
end
nBC = numel(PDE.BC);

% Which variables does each component depend on, and what are the orders
x_vars = logical(PDE.x_tab(:,3:2+N));
w_vars = logical(PDE.w_tab(:,3:2+N));
u_vars = logical(PDE.u_tab(:,3:2+N));
x_diff = PDE.x_tab(:,2+N+(1:N));
dep_vars = struct('x',x_vars,'w',w_vars,'u',u_vars);
var_names = cell(N,2);
for j=1:N
    var_names{j,1} = PDE.vars(j,1).varname{1};
    var_names{j,2} = PDE.vars(j,2).varname{1};
end


% % % Determine to which direction each condition belongs, and check that
% % % it is a condition that can be imposed along this direction (NOTES)
types = {'x','w','u'};
BC_dir = zeros(nBC,1);
for ii=1:nBC
    BC = PDE.BC{ii};
    BC_name = ['BC{',num2str(ii),'}'];
    varies_in = logical(PDE.BC_tab(ii,3:2+N));      % Variables in which the BC function varies

    % % Determine which variables are fixed by the terms in the condition
    fixed = false(1,N);
    for jj=1:numel(BC.term)
        term = BC.term{jj};
        obj_type = types{find(isfield(term,types),1)};
        comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
        pos = cumsum(comp_vars);
        for j=find(comp_vars)
            kk = pos(j);
            % A variable is fixed if the term is evaluated at a boundary,
            % or integrated over the full domain.
            if strcmp(obj_type,'x') && (isa(term.loc(kk),'double') || isdouble(polynomial(term.loc(kk))))
                fixed(j) = true;
            elseif numel(term.I)>=kk && ~isempty(term.I{kk})
                I_lim = polynomial(term.I{kk});
                if isdouble(I_lim)
                    I_lim = double(I_lim);
                    fixed(j) = fixed(j) || (numel(I_lim)==2 && isequal(sort(I_lim(:))',sort(PDE.dom(j,:))));
                end
            end
        end
    end
    if ~any(fixed)
        error('split_BVPs:noBoundary',['The condition ',BC_name,' does not evaluate any state at a boundary, and is not supported as a boundary condition along any spatial direction.'])
    elseif nnz(fixed)>1
        error('split_BVPs:multipleDirections','Boundary conditions that fix several variables at once are only supported if every term of a state with variables does so, without inputs.')
    end
    i = find(fixed);
    if varies_in(i)
        error('split_BVPs:mixedDirection',['The condition ',BC_name,' involves a term that is evaluated at a boundary in ',var_names{i,1},', but also a term that varies in ',var_names{i,1},'; it cannot be separated into conditions along a single spatial direction.'])
    end

    % % Check the terms
    % Determine which variables the states that depend on s_i vary in
    vary_vars = [];
    for jj=1:numel(BC.term)
        term = BC.term{jj};
        obj_type = types{find(isfield(term,types),1)};
        comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
        if strcmp(obj_type,'x') && comp_vars(i)
            if isfield(term,'delay') && any(term.delay(:)~=0)
                error('split_BVPs:delay',['The condition ',BC_name,' involves a delay of a state; delays are not supported.'])
            end
            if isempty(vary_vars)
                vary_vars = comp_vars;
            elseif ~isequal(vary_vars,comp_vars)
                error('split_BVPs:mixedDomains',['The condition ',BC_name,' couples states that depend on ',var_names{i,1},' but vary in different sets of variables; this would couple the problems along ',var_names{i,1},' for different values of the other variables.'])
            end
        end
    end
    if isempty(vary_vars)
        error('split_BVPs:noState',['The condition ',BC_name,' does not involve any state that depends on ',var_names{i,1},', and cannot be a condition along this direction.'])
    end
    if any(varies_in & ~vary_vars)
        error('split_BVPs:extraVars',['The condition ',BC_name,' varies in variables ',strjoin(var_names(varies_in & ~vary_vars,1)',', '),', on which the states it constrains do not depend.'])
    end

    % Check the derivatives and coefficients of the states depending on s_i
    % (all other terms are moved to the forcing, which is not restricted)
    for jj=1:numel(BC.term)
        term = BC.term{jj};
        obj_type = types{find(isfield(term,types),1)};
        comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
        pos = cumsum(comp_vars);
        if strcmp(obj_type,'x') && comp_vars(i)
            % Derivatives of boundary values along other directions are not
            % separable
            D_tangent = term.D(pos(comp_vars & ~((1:N)==i)));
            if any(D_tangent)
                error('split_BVPs:tangentialDerivative',['The condition ',BC_name,' involves a derivative of a state along a direction other than ',var_names{i,1},'; this is not supported.'])
            end
            % The coefficients cannot vary along the boundary, as they
            % define the 1D problem, except on the variable of integration
            % of an integral over s_i (s_i or its dummy variable)
            allowed = cell(1,0);
            if numel(term.I)>=pos(i) && ~isempty(term.I{pos(i)})
                allowed = var_names(i,:);
            end
            if isa(term.C,'polynomial') && ~isdouble(term.C) && any(~ismember(term.C.varname,allowed))
                error('split_BVPs:varyingCoefficient',['The condition ',BC_name,' involves a state depending on ',var_names{i,1},' with coefficients that vary along the boundary; this is not supported.'])
            end
        end
    end
    BC_dir(ii) = i;
end


% % % Write the conditions along each direction s_i in terms of the states
% % % depending on s_i, regrouped by their order in s_i (u=[u_0;...;u_d]),
% % % as E*[D_{j}u_{j}(a);D_{j}u_{j}(b)]_{j} + int F(s)*u(s) ds = g, where
% % % the right-hand side g is the sum of the other terms (the forcing)
dir_error = @(i,msg) error('convert_PIETOOLS_BVP:direction','%s',['Boundary conditions along ',var_names{i,1},': ',msg]);
dirs = struct('BC_idcs',cell(1,N),'nrow',[],'x_idcs',[],'n',[],'x_order',[],'E',[],'F',[],'lab_E',[],'lab_F',[],'forcing',[]);
for i=1:N
    BC_idcs = find(BC_dir==i)';
    nrow = reshape(cellfun(@(c) c.size,PDE.BC(BC_idcs)),1,[]);
    x_idcs = find(dep_vars.x(:,i))';
    x_new_idx = zeros(numel(PDE.x),1);      x_new_idx(x_idcs) = 1:numel(x_idcs);
    nk = PDE.x_tab(x_idcs,2);   dk = x_diff(x_idcs,i);      % sizes and orders
    d = max([0;dk]);
    dirs(i).BC_idcs = BC_idcs;      dirs(i).nrow = nrow;      dirs(i).x_idcs = x_idcs;
    dirs(i).forcing = cell(1,numel(BC_idcs));
    if d==0
        if ~isempty(BC_idcs)
            dir_error(i,'None of the states is differentiable, so that there are no boundary values to be constrained.')
        end
        continue
    end
    % Regroup the states by order, keeping the order within each group;
    % pos(k) is the offset of state k in u_{dk(k)}
    [~,x_order] = sort(dk(:)');
    n = zeros(1,d+1);       pos = zeros(numel(x_idcs),1);
    for k=1:numel(x_idcs)
        pos(k) = n(dk(k)+1);
        n(dk(k)+1) = n(dk(k)+1)+nk(k);
    end
    n_tot = sum((1:d).*n(2:end));               % number of boundary values
    row0 = [0,cumsum(nrow)];
    if row0(end)<n_tot
        dir_error(i,['The number of boundary conditions (',num2str(row0(end)),') should be at least the total order of the states (',num2str(n_tot),').'])
    end
    % Columns of E: [D_{j}u_{j}(a);D_{j}u_{j}(b)] for each order j, of F:
    % u_{j}; lab_E, lab_F give the state of each column
    idx_E0 = cumsum([0,2*(1:d).*n(2:end)]);     % start of block j in E
    idx_F0 = cumsum([0,n]);                     % start of block j in F
    E = zeros(row0(end),2*n_tot);       F = polynomial(zeros(row0(end),sum(n)));
    lab_E = zeros(1,2*n_tot);           lab_F = zeros(1,sum(n));
    for k=1:numel(x_idcs)
        j = dk(k);      idx = pos(k)+(1:nk(k));
        lab_F(idx_F0(j+1)+idx) = k;
        for dd=0:j-1
            lab_E(idx_E0(j)+dd*n(j+1)+idx) = k;
            lab_E(idx_E0(j)+(j+dd)*n(j+1)+idx) = k;
        end
    end
    for kk=1:numel(BC_idcs)
        rws = row0(kk)+(1:nrow(kk));
        Rsel = zeros(row0(end),nrow(kk));   Rsel(rws,:) = eye(nrow(kk));
        dirs(i).forcing{kk} = cell(1,0);
        for jj=1:numel(PDE.BC{BC_idcs(kk)}.term)
            term = PDE.BC{BC_idcs(kk)}.term{jj};
            obj_type = types{find(isfield(term,types),1)};
            comp_vars = dep_vars.(obj_type)(term.(obj_type),:);
            if ~strcmp(obj_type,'x') || ~comp_vars(i)
                % Any other term goes to the right-hand side
                term.C = -term.C;
                dirs(i).forcing{kk}{end+1} = term;
                continue
            end
            % A state depending on s_i: its boundary value (derivative
            % D<j at a or b) enters E, its integral over [a,b] enters F
            ki = nnz(comp_vars(1:i));               % position of s_i in its variables
            k = x_new_idx(term.x);      j = dk(k);      cols = pos(k)+(1:nk(k));
            C = polynomial(term.C);
            if isempty(term.I{ki})
                isb = double(polynomial(term.loc(ki)))==PDE.dom(i,2);
                c = idx_E0(j) + isb*j*n(j+1) + term.D(ki)*n(j+1) + cols;
                E(rws,c) = E(rws,c) + double(C);
            else
                if ismember(var_names{i,2},C.varname)
                    C = subs(C,PDE.vars(i,2),PDE.vars(i,1));
                end
                % An integral of a derivative D^{r}x is integrated by parts
                % into boundary values (E) and an integral of x (F):
                %   int C*D^{r}x = sum_{m<r} (-1)^m [C^{(m)}*D^{r-1-m}x]_a^b
                %                  + (-1)^r int C^{(r)}*x
                r = term.D(ki);
                for m=0:r-1                                                 % DJ (CC), 11/13/2026
                    c = idx_E0(j) + (r-1-m)*n(j+1) + cols;      % D^{r-1-m}x(a)
                    E(rws,c) = E(rws,c) - (-1)^m*double(subs(C,PDE.vars(i,1),PDE.dom(i,1)));
                    E(rws,c+j*n(j+1)) = E(rws,c+j*n(j+1)) + (-1)^m*double(subs(C,PDE.vars(i,1),PDE.dom(i,2)));
                    C = diff(C,PDE.vars(i,1));
                end
                Csel = zeros(nk(k),sum(n));     Csel(:,idx_F0(j+1)+cols) = eye(nk(k));
                F = F + Rsel*(-1)^r*C*Csel;
            end
        end
    end
    dirs(i).n = n;      dirs(i).x_order = x_order;      dirs(i).E = E;      dirs(i).F = F;
    dirs(i).lab_E = lab_E;      dirs(i).lab_F = lab_F;
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
