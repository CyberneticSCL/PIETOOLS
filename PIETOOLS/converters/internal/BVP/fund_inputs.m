function [inputs_ops,inputs_info] = fund_inputs(group,dirs)
% [inputs_ops,inputs_info] = FUND_INPUTS(group,dirs) parametrizes
% the inputs in the boundary conditions of a group of states by free
% functions wf_c ("fundamental inputs", uf_c for actuator inputs),
% such that every wf_c gives a state in the PDE domain, and every state
% in the domain with compatible boundary data is obtained from some
% wf_c.
%
% INPUT
% - group:  struct describing the q states x_pde_k, with ns scalar
%           components in total, that depend on the variables s_1,...,s_N,
%           with fields
%       vars:       1 x N cellstr, the names of the variables;
%       dom:        N x 2 array, the interval of each variable;
%       ns:         the number of scalar components of x_pde;
%       Delta:      q x N array, the order up to which each state is
%                   differentiable in each variable;
%       alphas:     nal x N array, the derivative orders alpha for which
%                   D^alpha x_pde is needed;
%       stpos:      ns x 1 array, the state of each component of x_pde;
%       cs_start:   1 x q array, the position in x_pde of the first
%                   component of each state;
% - dirs:   1 x N struct array, describing the 1D boundary value problem
%           along each direction s_i, with fields
%       triv:       true if there are no conditions along s_i, or if s_i is
%                   passed through (rank deficient, for inputs derived from
%                   other states, see compose_Tmaps);
%       forced:     true if some conditions along s_i have inputs;
%       R1:         cell, with R1{a+1} the 'sopvar' mapping the 1D solution
%                   to the a-th derivative of x_pde in s_i (R1{1} is
%                   T_{1,i});
%       Rp:         cell, with Rp{a+1} the same for the Taylor polynomial
%                   (in s_i) of the boundary data c;
%       Cw:         matrix with c=Cw*g, for c the boundary data and g the
%                   right-hand sides of the conditions along s_i;
%       rows_in:    logical, the conditions (rows of g) that have inputs;
%       d,pos,ell:  the entry of c in position k is the Taylor coefficient
%                   ell(k) (of order 0,...,d-1) at the lower boundary of
%                   the component pos(k) of x_pde;
%       src:        struct with fields pde, elem, scale, dep (arrays with a
%                   row for each condition): the condition has the
%                   right-hand side scale*w_{pde}(elem), a single element
%                   of a single input that depends on the variables dep
%                   (logical, N columns) among those other than s_i;
%                   pde=0 for other conditions;
%
% OUTPUT
% - inputs_ops:   nW x nal cell, with inputs_ops{c,a} the operator
%           wf_c --> D^{alphas(a,:)}x_pde;
% - inputs_info:  1 x nW struct array describing each wf_c, see
%           info.inputs of compose_Tmaps;
%
% EXAMPLE
% For one state u(s1,s2) on [0,1]^2, second order in both variables, with
% u(0,s2)=w1(s2), u(s1,0)=w2(s1), u(1,s2)=u(s1,1)=0, there are three
% fundamental inputs: wf_1=w1(0), a number (the corner s1=s2=0),
% wf_2(s1)=w2''(s1) and wf_3(s2)=w1''(s2). The right-hand sides w1
% and w2 are combinations of these, so that w1(0)=w2(0) holds
% automatically. With also u(1,s2)=w3(s2) and u(s1,1)=w4(s1) inputs, the
% corner data are the four numbers w1(0), w1'(0), w3(0) and w3'(0): the
% Taylor coefficients at the lower boundary of the inputs along s2 on
% s1=0 and s1=1. The other inputs are the second derivatives of w1,...,w4.
%
% NOTES
% - The boundary data along a direction with inputs are the Taylor
%   coefficients of the state at its lower boundary. Data along different
%   directions are tied together (e.g. at corners), so they are not all
%   free. Instead, ghat_J is defined for each set J of the directions with
%   conditions that contains those without inputs, and not all of them. It
%   is a function of the variables in J (and of those without conditions),
%   with values the data along the other directions, restricted to those
%   allowed by the conditions without inputs.
% - The operators apply, from the last direction to the first, the
%   derivative of T_{1,j} (j in J), or of the Taylor polynomial of the data
%   (other directions).
% - An input that does not depend on all variables of its boundary is
%   constant in the others. The pieces with a derivative in such a
%   variable vanish (forced), and the others are functions of the variables
%   of the component that the input depends on, extended as constants (E
%   in the main function). The pieces that vanish are chosen first, so the
%   others are independent of them, and the right-hand sides are tested to
%   be constant at the end, for all fundamental inputs. E.g. for scalar
%   inputs on all four faces of a 2D state, the corner conditions make
%   them equal, and there is one fundamental input.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - fund_inputs
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
% DJ (CC), 11/13/2026: Drop the output rhs_ops;

% Directions: Td without conditions, Kd with conditions, Fd with inputs,
% Hd with conditions without inputs (homogeneous, preserved by the T_{1,j})
N = numel(group.vars);      ns = group.ns;      nal = size(group.alphas,1);
Td = find([dirs.triv]);     Kd = find(~[dirs.triv]);
Fd = find([dirs.forced]);   Hd = Kd(~[dirs(Kd).forced]);
nF = numel(Fd);

% % % Maps relating the boundary data c along s_i (i in Fd) to the
% conditions c=Cw*g, added to dirs(i) as the fields
%   cidx:   d x ns array, the position in c of each (Taylor coefficient,
%           component of x_pde) pair;
%   Cinv:   inverse of Cw, so that g=Cinv*c;
%   Gbox:   maps the box R^d x R^ns of Taylor coefficients to the
%           right-hand sides g of the conditions with inputs;
%   Sj:     cell, Sj{a+1} the a-th derivative of the Taylor polynomial of
%           the data, acting on the box;
%   Mu:     matrix such that Mu*c=0 must hold: the rows of Cinv of the
%           conditions without inputs, with redundant rows removed;
for i=Fd
    ntot = numel(dirs(i).pos);
    Cw_inv = inv(dirs(i).Cw);
    rows_inp = dirs(i).rows_in;     rows_hom = find(~dirs(i).rows_in);
    cidx = zeros(dirs(i).d,ns);
    cidx(sub2ind([dirs(i).d,ns],dirs(i).ell+1,dirs(i).pos)) = 1:ntot;
    embed_box = sparse(1:ntot,dirs(i).ell*ns+dirs(i).pos,1,ntot,dirs(i).d*ns);
    dirs(i).cidx = cidx;
    dirs(i).Cinv = Cw_inv;
    dirs(i).Gbox = Cw_inv(rows_inp,:)*embed_box;
    dirs(i).Sj = cell(size(dirs(i).Rp));
    for a=find(~cellfun(@isempty,dirs(i).Rp(:)))'
        dirs(i).Sj{a} = dirs(i).Rp{a}*embed_box;
    end
    % The conditions without inputs must act on one class of states (same
    % orders in the other directions) each, so that they commute with the
    % Taylor polynomials; redundant rows are removed per class (SVD)
    Mu = Cw_inv(rows_hom,:);
    if ~isempty(rows_hom)
        oth_idcs = [1:i-1,i+1:N];
        if isempty(oth_idcs)
            class_id = ones(size(group.Delta,1),1);
        else
            [~,~,class_id] = unique(group.Delta(:,oth_idcs),'rows');
        end
        class_pos = class_id(group.stpos(dirs(i).pos));    class_pos = class_pos(:)';
        tol = 1e-10*max(1,norm(Mu,'fro'));
        Mu_new = zeros(0,ntot);
        for class_k=unique(class_pos)
            cols = class_pos==class_k;
            [~,~,V] = svd(Mu(:,cols),'econ');
            r = sum(svd(Mu(:,cols))>tol);
            Mu_class = zeros(r,ntot);       Mu_class(:,cols) = V(:,1:r)';
            Mu_new = [Mu_new;Mu_class];                                                  %#ok<AGROW>
        end
        if size(Mu_new,1)~=rank(Mu,tol)
            error('compose_Tmaps:inputRowsMixed',['The boundary conditions along ',group.vars{i},' without inputs couple components that have different orders in another direction, so that the inputs cannot be parameterized.'])
        end
        Mu = Mu_new;
    end
    dirs(i).Mu = Mu;
end

% One component ghat_J for each set J of Kd containing Hd, but not all of
% Kd (s enumerates the subsets of Fd). The boundary data along Ic=Fd\J are
% the values of ghat_J, and I=J and Td its variables.
% The values of ghat_J are the data (normalized conditions) along Ic, in
% free coordinates Z; blocks lists these components with their operators
blocks = struct('J',{},'I',{},'Ic',{},'dim',{},'ops',{},'rhs',{},'parts',{},'B',{});
% The lifted 1D operators recur in many components and derivatives, so they
% are kept in cache
cache = containers.Map('KeyType','char','ValueType','any');
for s=0:2^nF-2
    inJ = logical(bitget(s,1:nF));
    J = sort([Hd,Fd(inJ)]);     Ic = Fd(~inJ);      I = sort([J,Td]);

    % % % The entries of the values of ghat_J, and a basis Z of the values
    % that the conditions without inputs allow
    % ell (mu x p): the Taylor coefficient of each entry in each of the p
    % directions of Ic, up to the order of the state; pos (mu x 1): the
    % component of x_pde of each entry. Z (mu x r): the vectors of values
    % with Mu*c=0 in each direction of Ic, for the entries (l_t,pos) with
    % the other coefficients fixed
    p = numel(Ic);      dIc = arrayfun(@(k) dirs(k).d,Ic);
    % All (l_1,...,l_p) in the box R^{d_1} x...x R^{d_p}, the last fastest
    ells = zeros(1,0);
    for t=1:p
        ells = [repmat(ells,dIc(t),1),kron((0:dIc(t)-1)',ones(size(ells,1),1))];
    end
    [ell_idx,pos_idx] = ndgrid(1:size(ells,1),1:ns);    ell_idx = ell_idx(:);     pos_idx = pos_idx(:);
    % Keep the entries (l_1,...,l_p,pos) with l_t < the order of the state
    % of pos in direction i_t
    valid = true(numel(ell_idx),1);
    for t=1:p
        valid = valid & ells(ell_idx,t)<group.Delta(group.stpos(pos_idx),Ic(t));
    end
    ell = ells(ell_idx(valid),:);        pos = pos_idx(valid);        mu = numel(pos);
    if mu==0
        continue
    end
    % Constraints: for each t, Mu of direction i_t applied to the entries
    % (l_t,pos), for every choice of the other coefficients
    con_keys = zeros(0,p+1);    con_cols = zeros(0,1);        con_vals = zeros(0,1);
    for t=1:p
        Mu = dirs(Ic(t)).Mu;
        if isempty(Mu)
            continue
        end
        cidx_ent = dirs(Ic(t)).cidx(sub2ind(size(dirs(Ic(t)).cidx),ell(:,t)+1,pos));
        [mu_row,entry_idx,mu_val] = find(Mu(:,cidx_ent));
        mu_row = mu_row(:);     entry_idx = entry_idx(:);     mu_val = mu_val(:);
        oth_ell = ell(:,[1:t-1,t+1:p]);
        con_keys = [con_keys;[t*ones(numel(mu_row),1),mu_row,oth_ell(entry_idx,:)]];                      %#ok<AGROW>
        con_cols = [con_cols;entry_idx];       con_vals = [con_vals;mu_val];                                        %#ok<AGROW>
    end
    % Z: null space of the stacked constraint matrix
    Z = eye(mu);
    if ~isempty(con_vals)
        [~,~,row_id] = unique(con_keys,'rows');
        Cmat = full(sparse(row_id,con_cols,con_vals,max(row_id),mu));
        [~,~,V] = svd(Cmat);
        sv = svd(Cmat);
        Z = V(:,sum(sv>1e-10*max(1,max(sv)))+1:end);
        Z(abs(Z)<1e-12) = 0;
    end
    if size(Z,2)==0
        continue
    end

    [ops_b,rhs_b] = component_operators(group,dirs,cache,J,I,Ic,ell,pos,Z);
    [parts,B] = choose_pieces(group,dirs,Ic,J,ell,pos,Z);
    blocks(end+1) = struct('J',J,'I',I,'Ic',Ic,'dim',size(Z,2),'ops',{ops_b},'rhs',{rhs_b},'parts',{parts},'B',B);   %#ok<AGROW>
end
nb = numel(blocks);

% % % Change to the coordinates that are derivatives of inputs of the PDE
% The new coordinates of block b are the pieces
%       p_b = B_b*y_b + sum_b2 L_{b,b2}*y_b2
% for the old coordinates y, with b2 the blocks whose set J strictly
% contains that of b (the other blocks do not contribute to the
% derivatives and boundary values of the right-hand sides, as the data
% along a direction of J_b that is not in J_b2 are polynomials of lower
% degree than the order of the derivative). Solving from the largest J
% down, y_b=sum_b3 X_{b,b3}*p_b3 for the new coordinates p_b3, and the
% operators on y are composed with X. rhs_ops{ii,c} maps wf_c to the rows
% with inputs of the BCs along the ii-th direction with inputs ([] if
% wf_c does not depend on it), for the check at the end.
inputs_info = struct('vars',{},'dom',{},'dim',{},'pde',{},'Dorder',{},'loc',{},'sel',{});
inputs_ops = cell(0,nal);       rhs_ops = cell(nF,0);
size_J = arrayfun(@(blk) numel(blk.J),blocks);
Lmaps = cell(nb,nb);
for b3=1:nb
    % The pieces that are not forced to vanish (see choose_pieces) are the
    % new coordinates, the others are zero
    forced = cell2mat(arrayfun(@(pt) repmat(pt.forced,numel(pt.elem),1),blocks(b3).parts(:),'UniformOutput',false));
    if all(forced)
        continue
    end
    X = cell(1,nb);     Binv = inv(blocks(b3).B);       X{b3} = Binv(:,~forced);
    inside = find(arrayfun(@(blk) all(ismember(blk.J,blocks(b3).J)),blocks) & size_J<size_J(b3));
    [~,order] = sort(size_J(inside),'descend');
    for b=inside(order)
        acc = [];
        for b2=find(~cellfun(@isempty,X))
            if size_J(b2)<=size_J(b) || ~all(ismember(blocks(b).J,blocks(b2).J))
                continue
            end
            if isempty(Lmaps{b,b2})
                % The map L from the old coordinates of b2 to the pieces
                % of b, a row for each piece (in the order of the parts):
                % the derivative of the right-hand side of the condition
                % of the piece, of the order of the piece, evaluated at
                % the lower boundary in the other directions of Ic, divided
                % by the coefficient of the input in it. The right-hand
                % side is the output of the operator rhs of b2, so that
                % its derivative and value are diff and subs of that
                % 'sopvar'. The right-hand sides along a direction in J of
                % b2 do not depend on b2 (zero).
                vars_in = group.vars(blocks(b2).I);     vars_out = group.vars(blocks(b).I);
                Lrows = cell(blocks(b).dim,1);          p = 0;
                for g=1:numel(blocks(b).parts)
                    part = blocks(b).parts(g);
                    for k=1:size(part.rows,1)
                        p = p+1;
                        i = part.rows(k,1);         rho = part.rows(k,2);       scale = part.rows(k,3);
                        G = blocks(b2).rhs{Fd==i};
                        if isempty(G)
                            Lrows{p} = mat2sopvar(zeros(1,blocks(b2).dim),struct('in',{vars_in},'out',{vars_out}),struct('in',group.dom(blocks(b2).I,:),'out',group.dom(blocks(b).I,:)));
                            continue
                        end
                        Row = G(find(find(dirs(i).rows_in)==rho),:);
                        dv = find(part.diff>0);
                        if ~isempty(dv)
                            Row = diff(Row,group.vars(dv),part.diff(dv));
                        end
                        ev = blocks(b).Ic(blocks(b).Ic~=i);
                        if ~isempty(ev)
                            Row = subs(Row,group.vars(ev),group.dom(ev,1)');
                        end
                        Lrows{p} = Row*(1/scale);
                    end
                end
                Lmaps{b,b2} = vertcat(Lrows{:});
            end
            if isempty(acc)
                acc = Lmaps{b,b2}*X{b2};
            else
                acc = acc+Lmaps{b,b2}*X{b2};
            end
        end
        if ~isempty(acc)
            X{b} = -inv(blocks(b).B)*acc;
        end
    end
    % The operators of each part (the new coordinates of the same input
    % and derivative), as the old operators composed with X
    first = 0;
    for g=find(~[blocks(b3).parts.forced])
        part = blocks(b3).parts(g);
        cols = first+(1:numel(part.elem));      first = first+numel(cols);
        c = size(inputs_ops,1)+1;
        bs = find(~cellfun(@isempty,X));
        % An input that does not depend on all variables of the component
        % is constant in the others: the new coordinates are functions of
        % the variables W of the component that the input depends on,
        % extended as constants by E
        I = blocks(b3).I;       W = I(part.dep(I));
        if numel(W)<numel(I)
            E = mat2sopvar(eye(numel(cols)),struct('in',{group.vars(W)},'out',{group.vars(I)}),struct('in',group.dom(W,:),'out',group.dom(I,:)));
        end
        for a=1:nal
            inputs_ops{c,a} = blocks(bs(1)).ops{a}*X{bs(1)}(:,cols);
            for b=bs(2:end)
                inputs_ops{c,a} = inputs_ops{c,a}+blocks(b).ops{a}*X{b}(:,cols);
            end
            if numel(W)<numel(I)
                inputs_ops{c,a} = inputs_ops{c,a}*E;
            end
        end
        % The right-hand sides: blocks without the operator add nothing
        for ii=1:nF
            rhs_ops{ii,c} = [];
            for b=bs
                if isempty(blocks(b).rhs{ii})
                    continue
                end
                term = blocks(b).rhs{ii}*X{b}(:,cols);
                if isempty(rhs_ops{ii,c})
                    rhs_ops{ii,c} = term;
                else
                    rhs_ops{ii,c} = rhs_ops{ii,c}+term;
                end
            end
            if ~isempty(rhs_ops{ii,c}) && numel(W)<numel(I)
                rhs_ops{ii,c} = rhs_ops{ii,c}*E;
            end
        end
        % Record the input: its variables, and the derivative of the input
        % of the PDE that it is
        inputs_info(c) = struct('vars',{group.vars(W)},'dom',group.dom(W,:),'dim',numel(cols),'pde',part.pde,'Dorder',part.diff,'loc',part.loc,'sel',part.elem(:));
    end
end

% % % Check that the right-hand sides of the conditions with an input that
% does not depend on all variables of the boundary are constant in the
% others, for all fundamental inputs. The pieces that vanish are set to
% zero, which gives all data compatible with the conditions if the other
% conditions that the constancy implies do not restrict the other pieces.
% If they do (e.g. a constant input on a boundary that also needs a
% derivative along it to equal that of another input), an error is thrown
for ii=1:nF
    i = Fd(ii);     rows_i = find(dirs(i).rows_in);
    for k=1:numel(rows_i)
        U = find(~dirs(i).src.dep(rows_i(k),:));        U = U(U~=i);
        for c=1:size(rhs_ops,2)
            if isempty(rhs_ops{ii,c}) || isempty(U)
                continue
            end
            for u=U
                d = diff(rhs_ops{ii,c}(k,:),group.vars(u),1);
                if ~all(cellfun(@(P) nnz(abs(P)>1e-10)==0,d.params))
                    error('compose_Tmaps:inputsNotRepresentable',['The inputs in the boundary conditions along ',group.vars{i},' that do not depend on ',group.vars{u},' are not constant in ',group.vars{u},' for all compatible data: the conditions along the other directions impose relations between the inputs that are not supported.'])
                end
            end
        end
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [parts,B] = choose_pieces(group,dirs,Ic,J,ell,pos,Z)
% [parts,B] = CHOOSE_PIECES(GROUP,DIRS,IC,J,ELL,POS,Z) chooses coordinates
% of the values of ghat_J that are derivatives of inputs of the PDE.
%
% INPUT
% - group,dirs:   structs, see fund_inputs, with the fields added there to
%                 dirs;
% - Ic,J:         the directions Fd\J of the data in ghat_J, and J;
% - ell,pos,Z:    the entries and the basis of their values, see
%                 fund_inputs;
%
% OUTPUT
% - parts:  struct array, one for each derivative of an input of the PDE
%           (the pieces with the same input, derivative and location are
%           one part), with fields
%           - pde:  index of the input of the PDE (in the group);
%           - diff: 1 x N array, the orders of its derivative;
%           - loc:  1 x N array, -1 for the variables of the input in which
%                   it is evaluated at the lower boundary, 0 otherwise;
%           - dep:  1 x N logical array, the variables of the input;
%           - forced:   true if the piece vanishes: the derivative is in a
%                   variable that the input does not depend on;
%           - elem: array, the element of the input of each piece;
%           - rows: array with a row [i,row,scale] for each piece, the
%                   direction and row of the condition with the input, and
%                   the coefficient of the element in it;
% - B:      dim x dim invertible matrix, giving the part of each piece (in
%           the order of the parts) that is in the values of ghat_J itself,
%           as B*(coordinates of the values in Z);
%
% NOTES
% - A piece is a derivative, at the lower boundary in the other directions
%   of Ic, of a single element of an input of the PDE (a candidate, see
%   below). The pieces are chosen greedily, as long as their parts in the
%   values are independent of those already chosen: the other parts of the
%   pieces, from the components with a larger J, are triangular (see
%   fund_inputs). The candidates that vanish (forced) are taken first, so
%   that the others are independent of them: the values of the conditions
%   that vanish then only restrict the pieces that are forced, which are
%   set to zero, and the others are free (if the vanishing conditions do
%   not restrict the components with a larger J, which is checked at the
%   end of fund_inputs).
% - The candidates are the derivatives of an element of an input of the
%   PDE that have a part in the values of ghat_J, with the coefficients
%   coef of this part in the free coordinates Z. Let g be the right-hand
%   side of a condition with an input along s_i (i in Ic). The data along
%   s_i are c=Cw*g, so g is Cinv(row,:)*c. The values of ghat_J are the
%   data along the other directions k of Ic, which are Taylor coefficients
%   (ell_k) at the lower boundary of the polynomial in s_k that the values
%   contribute to c, so the part of the Taylor coefficient lp (and the
%   derivative of the states in the variables of J) of g in the values is
%   sum_{ell_i} Cinv(row,ell_i)*value(ell_i,lp).
%   A candidate has fields
%   - pde,elem: the input of the PDE, and its element;
%   - dir,row,scale: the direction i in Ic and row of the condition along
%               s_i with the input, with right-hand side
%               scale*w_{pde}(elem);
%   - diff:     1 x N array, the orders of the derivative of the input: the
%               order of the states in the variables of J, and the Taylor
%               order in the other variables of Ic;
%   - loc:      1 x N array, -1 for the variables of Ic other than i that
%               the input depends on, in which it is evaluated at the lower
%               boundary, 0 otherwise;
%   - dep:      1 x N logical array, the variables of the input;
%   - forced:   true if the derivative is in a variable that the input does
%               not depend on, so that it vanishes;
%   - coef:     1 x dim array, the part of the candidate in the values of
%               ghat_J, coef*(coordinates in Z);
N = numel(group.vars);      p = numel(Ic);      dim = size(Z,2);
cand = struct('pde',{},'elem',{},'dir',{},'row',{},'scale',{},'diff',{},'loc',{},'dep',{},'forced',{},'coef',{});
for t=1:p
    i = Ic(t);      oth = [1:t-1,t+1:p];
    % Column of c_i of each entry, and the Taylor orders in the others
    cols = dirs(i).cidx(sub2ind(size(dirs(i).cidx),ell(:,t)+1,pos));
    ell_oth = ell(:,oth);
    lps = zeros(1,0);
    if ~isempty(oth)
        lps = unique(ell_oth,'rows');
    end
    for rho=find(dirs(i).rows_in(:) & dirs(i).src.pde(:)>0)'
        for il=1:size(lps,1)
            lp = lps(il,:);
            a = dirs(i).Cinv(rho,cols)'.*all(ell_oth==lp,2);
            if ~any(abs(a)>1e-12)
                continue
            end
            % The derivative in the variables of J is that of the states
            diff_c = zeros(1,N);        loc_c = zeros(1,N);
            Dst = group.Delta(group.stpos(unique(pos(abs(a)>1e-12))),J);
            if ~isempty(J)
                if size(unique(Dst,'rows'),1)>1
                    continue
                end
                diff_c(J) = Dst(1,:);
            end
            diff_c(Ic(oth)) = lp;       loc_c(Ic(oth)) = -1;
            % The right-hand side is constant in the variables U of the
            % boundary that the input does not depend on
            dep = dirs(i).src.dep(rho,:);       U = true(1,N);      U(i) = false;       U(dep) = false;
            loc_c(U) = 0;
            scale = dirs(i).src.scale(rho);
            coef = (a'/scale)*Z;
            if ~any(abs(coef)>1e-12)
                continue
            end
            cand(end+1) = struct('pde',dirs(i).src.pde(rho),'elem',dirs(i).src.elem(rho),'dir',i,'row',rho,'scale',scale,'diff',diff_c,'loc',loc_c,'dep',dep,'forced',any(diff_c(U)>0),'coef',coef);   %#ok<AGROW>
        end
    end
end

chosen = zeros(1,0);        B = zeros(0,dim);
for k=[find([cand.forced]),find(~[cand.forced])]
    if rank([B;cand(k).coef],1e-10)>size(B,1)
        chosen(end+1) = k;          B(end+1,:) = cand(k).coef;                          %#ok<AGROW>
        if numel(chosen)==dim
            break
        end
    end
end
if numel(chosen)<dim
    error('compose_Tmaps:inputsNotRepresentable',['The inputs in the boundary conditions along ',strjoin(group.vars(Ic),', '),' cannot all be written as derivatives of inputs of the PDE: this needs each condition with an input to have a single element of a single input as right-hand side, depending on some of the variables other than the one the condition is along.'])
end
% Parts: the pieces with the same input, derivative and location
key = cell2mat(arrayfun(@(k) [cand(k).pde,cand(k).diff,cand(k).loc],chosen(:),'UniformOutput',false));
[~,~,part_id] = unique(key,'rows','stable');
[part_id,order] = sort(part_id);        chosen = chosen(order);         B = B(order,:);
parts = struct('pde',cell(1,max(part_id)),'diff',[],'loc',[],'dep',[],'forced',[],'elem',[],'rows',[]);
for g=1:max(part_id)
    in_part = chosen(part_id==g);
    rows = cell2mat(arrayfun(@(k) [cand(k).dir,cand(k).row,cand(k).scale],in_part(:),'UniformOutput',false));
    parts(g) = struct('pde',cand(in_part(1)).pde,'diff',cand(in_part(1)).diff,'loc',cand(in_part(1)).loc,'dep',cand(in_part(1)).dep,'forced',cand(in_part(1)).forced,'elem',[cand(in_part).elem],'rows',rows);
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [ghat_ops_row,rhs_ops_col] = component_operators(group,dirs,cache,J,I,Ic,ell,pos,Z)
% [ghat_ops_row,rhs_ops_col] = COMPONENT_OPERATORS(group,dirs,cache,J,I,Ic,
% ell,pos,Z) builds the operators of one component ghat_J: to the
% derivatives of the state, and to the right-hand sides of the conditions
% with inputs.
%
% INPUT
% - group,dirs:   structs, see fund_inputs, with the fields added there to
%                 dirs;
% - cache:        containers.Map in which the lifted 1D operators are
%                 kept, updated by this function;
% - J,I,Ic:       the directions of the set J, the variables I=J and Td
%                 of ghat_J, and the directions Ic=Fd\J of the data;
% - ell,pos,Z:    the entries and the basis of their values, see
%                 fund_inputs;
%
% OUTPUT
% - ghat_ops_row:  1 x nal cell, with ghat_ops_row{a} the operator
%           ghat_J --> D^{alphas(a,:)}x_pde;
% - rhs_ops_col:    nF x 1 cell, with rhs_ops_col{ii} the operator ghat_J
%           --> the right-hand sides with inputs of the conditions along
%           the ii-th direction with inputs, if that direction is in Ic,
%           and [] if it is in J;
vars = group.vars;          N = numel(vars);        ns = group.ns;
nal = size(group.alphas,1);
Fd = find([dirs.forced]);   nF = numel(Fd);
p = numel(Ic);              dIc = arrayfun(@(i) dirs(i).d,Ic);
mu = numel(pos);            boxsize = prod(dIc)*ns;
% Position in the box R^{d_1} x ... x R^{d_p} x R^{ns} (first factor
% slowest) of the entries, for the factors of the box having sizes box_dims
boxidx = @(ell,box_dims) (ell*[cumprod(box_dims(2:end),'reverse'),1]')*ns + pos;
Zbox = sparse(boxidx(ell,dIc),(1:mu)',1,boxsize,mu)*Z;

% Jobs 1,...,nal: ghat_J --> D^{alpha}x_pde. Jobs nal+ii: ghat_J --> g_i,
% the right-hand sides along s_i (i in Ic). Each job composes, from
% direction N to 1, the derivative of T_{1,j} (j in J) or of the Taylor
% polynomial (j in Ic).
ghat_ops_row = cell(1,nal);       rhs_ops_col = cell(nF,1);
% Job job: the map to D^{alpha}x_pde (job<=nal), or to a BC right-hand side
for job=1:nal+nF
    if job<=nal
        alpha_job = group.alphas(job,:);    box_dims = dIc;       skip = 0;       Zjob = Zbox;
    else
        ii = job-nal;                i = Fd(ii);
        if any(J==i)
            continue
        end
        t = find(Ic==i);            perm = [t,1:t-1,t+1:p];
        alpha_job = zeros(1,N);            box_dims = dIc(perm);         skip = i;
        Zjob = sparse(boxidx(ell(:,perm),box_dims),(1:mu)',1,boxsize,mu)*Z;
    end
    stages = struct('op',{},'var',{});
    nfactors = p;                                 % factors left in the box
    for j=N:-1:1
        if dirs(j).triv || j==skip
            continue
        end
        % The 1D operator of direction j, lifted to act on the first factor
        % left in the box (nfactors factors, of sizes box_dims(1:nfactors))
        % T_{1,j} for j in J, and the Taylor polynomial otherwise, which
        % consumes a factor
        if any(J==j)
            op_1D = dirs(j).R1{alpha_job(j)+1};        tag = 'R';
        else
            nfactors = nfactors-1;
            op_1D = dirs(j).Sj{alpha_job(j)+1};        tag = 'S';
        end
        key = sprintf('%s%d_%d_%d',tag,j,alpha_job(j),prod(box_dims(1:nfactors)));
        if ~isKey(cache,key)
            cache(key) = kron(speye(prod(box_dims(1:nfactors))),op_1D);
        end
        lifted_op = cache(key);
        stages(end+1) = struct('op',lifted_op,'var',vars{j});                          %#ok<AGROW>
    end
    job_op = compose_sopvar(stages,vars(I),vars,group.dom,boxsize)*Zjob;
    if job<=nal
        ghat_ops_row{job} = job_op;
    else
        rhs_ops_col{ii} = dirs(i).Gbox*job_op;
    end
end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
