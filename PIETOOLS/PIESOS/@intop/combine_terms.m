function Kop = combine_terms(Kop,degs,dims,depmat)
% KOP = COMBINE_TERMS(KOP,DEGS,DIMS,DEPMAT) combines the terms in the
% functional KOP that, when acting on the distributed monomial defined by
% DEGS, are defined in terms of the same integral.
%
% INPUTS
% - Kop:    'intop' object representing a functional in N spatial
%           variables, acting on a distributed monomial in which
%           Kop.depmat.in(j) factors depend on spatial variable j;
% - degs:   (optional) 1 x p array of nonnegative integers specifying the
%           degree of the p state variables in the distributed monomial, so
%           that if degs = [d1,...,dp], then
%               X = x1^{d1} o ... o xp^{dp}
%           where o is the tensor product and the power is a tensor product
%           power;
% - dims:   (optional) p x 1 array indicating the vector size of each of
%           the state variables;
% - depmat: (optional) p x N logical array, with element (i,j) indicating
%           whether state variable i depends on spatial variable j;
%
% OUTPUTS
% - Kop:    'intop' object representing the same functional as the input
%           when acting upon the specified distributed monomial;
%
% NOTES
% - If only Kop is specified, only the terms involving the same integral
%   are combined. Otherwise, all of degs, dims, and depmat must be
%   specified.
% - If depmat has more than N columns, the columns of the spatial variables
%   on which none of the states in the monomial depends are discarded.
% - Equal functionals are returned with equal intmat and kernels, so that
%   their coefficients can be compared cell by cell, and a functional that
%   is already combined is returned unchanged (idempotent).
% - Every integral is first written as a sum of cells, one for each
%   ordering of the dummy variables consistent with its limits. Each cell
%   is in canonical form: every variable is integrated from a (0) up to
%   the next larger variable, the largest over [a,b], and multipliers
%   [t,t] follow a variable that is not a multiplier. E.g. [0,3,0,3]
%   becomes [0,2,0,3] + [0,3,0,1], and [2,3,0,3] becomes [0,3,0,1]. The
%   rows of an 'intop' are ordered (see intop.canonicalize_intmat), so
%   that the limits are never reversed where a cell is nonzero, and the
%   signed integral and the integral of the indicator function agree.
%   Integrals over more than 8 variables are not expanded (a warning is
%   issued).
% - In each term, the copies of the same state are ordered by decreasing
%   order of integration (the last variable is integrated first, see
%   'int_order' in 'intop'), along the first spatial variable on which the
%   state depends, ties by the next variable, remaining ties keeping their
%   order. Copies of a state are interchangeable, so the terms that
%   coincide after this relabeling are added: [0,3,0,3] with kernel k
%   becomes [0,3,0,1] with kernel k(s1,s2)+k(s2,s1).
% - The kernel of a multiplier is written in terms of the variable at
%   which it is evaluated: a variable and the multipliers at it are
%   represented by the lowest of them.
% - Copies of a vector-valued state that are equal in all spatial
%   variables commute, x_i(s)x_j(s) = x_j(s)x_i(s), so only the sum of the
%   kernels on the columns (i,j) and (j,i) is determined. This sum is
%   placed on the column with nondecreasing indices.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - combine_terms
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
% DJ, 10/22/2026: Update comment on reversed limits, which do not occur;

if nargin==1
    % Just combine terms involving the same integral
    Kop = unique_intmat(Kop);
    return
end

% Make sure the operator is properly specified
if ~Kop.chkval && ~isvalid(Kop)
    error("The functional is not properly declared.")
end

% Check that the state indices are properly specified
if ~isnumeric(degs) || any(round(abs(degs))~=degs)
    error("Degrees should be specified as 1 x p array of nonnegative integers.")
end
% Establish the order of the independent variables in each term
degs = reshape(degs,1,[]);
p = numel(degs);        d = sum(degs);
state_idcs = repelem(1:p,degs);

% Check that the dimensions are properly specified
if ~isnumeric(dims) || numel(dims)~=p
    error("A vector dimension must be specified for each of the factors in the distributed monomials.")
end
dims = reshape(dims,1,[]);
dims_full = dims(state_idcs);
if prod(dims_full)~=size(Kop,2)
    error("Vector dimension of the distributed monomial does not match column dimension of the functional.")
end

% Check that the dependency array is appropriate
if size(depmat,2)~=numel(Kop.intmat)
    % Get rid of variables on which the state does not depend
    depmat = depmat(:,any(depmat(unique(state_idcs),:),1));
end
if size(depmat,2)~=numel(Kop.intmat) || size(depmat,1)~=p || ~islogical(depmat)
    error("Dependency array must be specified as p x N logical array, for p state variables and N spatial variables.")
end
depmat_in = depmat(state_idcs,:);
if any(sum(depmat_in,1)~=Kop.depmat.in)
    error("Input degrees of variables in functional do not match specified dependency array.")
end

% For each of the spatial variables, determine what index (degree) of dummy
% variable associated with this spatial variable appears in each of the
% factors in the distributed monomial
N = numel(Kop.intmat);
pvar_deg = zeros(d,N);
for i=1:N
    pvar_deg(depmat_in(:,i),i) = 1:nnz(depmat_in(:,i));
end


%% Step 1: Express every integral as a sum of cells
% One cell for each ordering of the dummy variables consistent with the
% limits, in canonical form (see NOTES). Equal functionals then have equal
% intmat, so that their coefficients can be compared cell by cell.
is_expd = false;
for i=1:N
    imat = Kop.intmat{i};
    if isempty(imat),   continue,   end         % no dummy variables
    % A row is a cell if the variables that are not multipliers form a
    % chain, each integrated from a up to the next, the largest over the
    % full domain, and each multiplier is evaluated at one of these
    is_ok = is_cell(imat);
    if all(is_ok)
        continue
    end
    % Expand the rows that are not cells; the others keep their parameters
    cells_i = num2cell(imat,2);
    sgn_i = num2cell(ones(size(imat,1),1));
    [urows,~,uidx] = unique(imat(~is_ok,:),'rows');
    ucells = cell(size(urows,1),1);     usgn = ucells;
    for k=1:size(urows,1)
        % Cells of the integral: one canonical row for each ordering of
        % the variables that are not multipliers in which each variable
        % lies between its limits (the diagonals have measure zero). The
        % kernel is the same on each cell, up to the sign: -1 if a
        % variable is integrated from the larger to the smaller limit.
        row = urows(k,:);
        nd = numel(row)/2;
        L = row(1:2:end);       U = row(2:2:end);
        is_dlt = L==U;
        % tgt: the variable at which each multiplier is evaluated
        tgt = 1:nd;
        for kk=1:nd
            tgt(is_dlt(tgt)) = L(tgt(is_dlt(tgt)));
        end
        fre = find(~is_dlt);        n = numel(fre);
        if n>8
            warning("Integral over %d variables is not expanded into cells.",n)
            ucells{k} = row;        usgn{k} = 1;
            continue
        end
        loc = zeros(1,nd);          loc(fre) = 1:n;
        % Orderings P(r,:) of the free variables, increasing, with
        % positions pos
        P = perms(1:n);
        [~,pos] = sort(P,2);
        ok_ord = true(size(P,1),1);
        sgn_ord = ones(size(P,1),1);
        for a=1:n
            lo = L(fre(a));         up = U(fre(a));
            if lo>=1 && up<=nd
                % Both limits are variables: the variable is between them,
                % in either order, and the sign is -1 if the lower limit is
                % the larger
                pl = pos(:,loc(tgt(lo)));       pu = pos(:,loc(tgt(up)));
                is_fwd = pl<pos(:,a) & pos(:,a)<pu;
                is_rev = pu<pos(:,a) & pos(:,a)<pl;
                ok_ord = ok_ord & (is_fwd | is_rev);
                sgn_ord(is_rev) = -1;
            elseif lo>=1
                % Only the lower limit is a variable (the upper is b)
                ok_ord = ok_ord & pos(:,loc(tgt(lo))) < pos(:,a);
            elseif up<=nd
                % Only the upper limit is a variable (the lower is a)
                ok_ord = ok_ord & pos(:,a) < pos(:,loc(tgt(up)));
            end
        end
        % ord(r,k): kth smallest variable in ordering r
        ord = fre(P(ok_ord,:));
        usgn{k} = sgn_ord(ok_ord);
        nr = size(ord,1);
        crows = zeros(nr,2*nd);
        for kk=1:n
            if kk<n
                up = ord(:,kk+1);
            else
                up = (nd+1)*ones(nr,1);
            end
            crows(sub2ind(size(crows),(1:nr)',2*ord(:,kk))) = up;
        end
        dl = find(is_dlt);
        crows(:,2*dl-1) = repmat(tgt(dl),nr,1);
        crows(:,2*dl) = repmat(tgt(dl),nr,1);
        ucells{k} = crows;
    end
    cells_i(~is_ok) = ucells(uidx);     sgn_i(~is_ok) = usgn(uidx);
    n_new = cellfun(@(a) size(a,1),cells_i);
    % New rows, each with the old row (parameters) it comes from
    src = reshape(repelem(1:size(imat,1),n_new'),[],1);
    rows_new = vertcat(cells_i{:});
    sgn_new = vertcat(sgn_i{:});
    if isempty(rows_new)
        continue
    end
    ip = cellfun(@(a) (1:size(a,1))',Kop.intmat,'UniformOutput',false);
    ip{i} = src;
    Kop.params = Kop.params(ip{:});
    % Integrating from the larger to the smaller limit changes the sign
    ineg = find(sgn_new<0);
    if ~isempty(ineg)
        ip{i} = ineg;
        pneg = Kop.params(ip{:});
        Kop.params(ip{:}) = cellfun(@(a) -a,pneg,'UniformOutput',false);
    end
    Kop.intmat{i} = rows_new;
    is_expd = true;
end
if is_expd
    % Add the parameters of cells that coincide
    Kop.chkval = true;
    Kop = unique_intmat(Kop);
end


%% Step 2: Order the copies of each state, relabeling the cells
% Loop over all state variables which appear with degree at least 2
for state_num = 1:p
    if degs(state_num)<2
        % We need at least two factors to switch their order
        continue
    end

    % Determine in what factors state k appears
    is_state_k = state_idcs==state_num;
    % Determine which spatial variables the state depends on
    deps_k = depmat(state_num,:);
    if ~any(deps_k)
        % State does not depend on any spatial variable -- no dummy
        % integration variables exist to reorder.
        continue
    end
    % Order the m copies of the state in each term by decreasing order of
    % integration (last variable first), along the first spatial variable
    % the state depends on, ties by the next variable, remaining ties
    % keeping their order. A term in this order is a fixed point, so that
    % combine_terms is idempotent.
    var_idcs = find(deps_k);    nv = numel(var_idcs);
    m = degs(state_num);
    int_order = Kop.int_order;
    % Classes of rows of intmat{var} with the same order of integration
    % of the m copies, as the dummy variables pvars_k of the variable
    int_lvl = cell(1,nv);       int_cls = cell(1,nv);
    for t=1:nv
        pvars_k = pvar_deg(is_state_k,var_idcs(t))';
        vt = var_idcs(t);
        % Order of integration of the variable each copy equals: itself,
        % or the variable a multiplier is evaluated at, so equal copies
        % tie (their order follows from the other spatial variables)
        imat = Kop.intmat{vt};
        is_dlt = imat(:,2*pvars_k-1)==imat(:,2*pvars_k);
        tgt_var = repmat(pvars_k,size(imat,1),1);
        lwr = imat(:,2*pvars_k-1);
        tgt_var(is_dlt) = lwr(is_dlt);
        io = int_order{vt};
        rr = repmat((1:size(io,1))',1,numel(pvars_k));
        int_order_t = io(sub2ind(size(io),rr,tgt_var));
        [int_lvl{t},~,int_cls{t}] = unique(int_order_t,'rows');
    end
    % Order of the copies for each combination of classes, i.e. for each
    % term. The terms are a product over the rows along each variable, so
    % empty terms (no nonzero parameters) cannot all be in order: they
    % keep theirs
    ncls = cellfun(@(a) size(a,1),int_lvl);
    swp_all = zeros(prod(ncls),m);
    mx = max(cellfun(@(a) max(a(:)),int_lvl));
    wts = ((mx+1).^(nv-1:-1:0)).';      % weights of lexicographic key
    is_nz = cellfun(@(a) nnz(a)>0,Kop.params);
    idx_nz = repmat({':'},1,N);
    sub = cell(1,nv);           key = zeros(m,nv);
    for c=1:prod(ncls)
        [sub{:}] = ind2sub([ncls,1],c);
        for t=1:nv
            key(:,t) = int_lvl{t}(sub{t},:)';
            idx_nz{var_idcs(t)} = find(int_cls{t}==sub{t});
        end
        if any(is_nz(idx_nz{:}),'all')
            % sort is stable: tied copies keep their order
            [~,swp_all(c,:)] = sort(key*wts,'descend');
        else
            swp_all(c,:) = 1:m;
        end
    end
    if all(swp_all==1:m,'all')
        % All terms are in order already: nothing to do
        continue
    end
    % Blocks of terms with the same order: products of sets of classes,
    % one set per spatial variable. Found by splitting the array gid of
    % orders along the first variable along which it varies, grouping
    % equal slices
    [swp_idcs,~,gid] = unique(swp_all,'rows');
    gid = reshape(gid,[ncls,1]);
    todo = {arrayfun(@(n) 1:n,ncls,'UniformOutput',false)};
    blks = cell(1,0);           blk_gid = zeros(1,0);
    while ~isempty(todo)
        cl = todo{end};         todo(end) = [];
        gsub = gid(cl{:});
        if all(gsub(:)==gsub(1))
            blks{end+1} = cl;               %#ok<AGROW>
            blk_gid(end+1) = gsub(1);       %#ok<AGROW>
            continue
        end
        for t=1:nv
            ord = [t,1:t-1,t+1:max(nv,2)];
            gsl = reshape(permute(gsub,ord),size(gsub,t),[]);
            [~,~,g] = unique(gsl,'rows');
            if max(g)>1
                break
            end
        end
        for j=1:max(g)
            cl_j = cl;          cl_j{t} = cl{t}(g==j);
            todo{end+1} = cl_j;             %#ok<AGROW>
        end
    end

    blk_ops = cell(1,numel(blks));
    tmp_idcs =[cellfun(@(a) 1:size(a,1),Kop.intmat,'UniformOutput',false),{1}];
    state_order_j = 1:numel(state_idcs);        % by default, keep same order of all variables
    state_order_k = state_order_j(is_state_k);
    for j=1:numel(blks)                 % one block of terms per order
        % Extract those terms of Kop for which we want to swap the order
        % of the states as per swp_idcs(blk_gid(j),:)
        Kop_tmp = Kop;
        for t=1:nv
            % Rows of the block along the variables the state depends on
            vt = var_idcs(t);
            rtn_idcs = find(ismember(int_cls{t},blks{j}{t}));
            Kop_tmp.intmat{vt} = Kop.intmat{vt}(rtn_idcs,:);
            tmp_idcs{vt} = rtn_idcs;
        end
        Kop_tmp.params = Kop.params(tmp_idcs{:});
        % Switch the order of the state variables to match the desired
        % order of the spatial variables
        swp_idcs_j = swp_idcs(blk_gid(j),:);
        state_order_j(is_state_k) = state_order_k(swp_idcs_j);
        blk_ops{j} = reorder_states(Kop_tmp,state_order_j,dims_full,depmat_in);
    end

    % Sum the blocks in one pass, instead of adding them one by one with
    % plus. The blocks derive from the same functional, so they share Zd,
    % vars and dom (Zd is not compared) and differ in intmat, params and the
    % order of the variables in Z.in.
    B = numel(blk_ops);
    if B==1
        Kop = blk_ops{1};
        continue
    end
    % Union of the monomials, in the order of first appearance
    Zu = cell(1,N);
    for i=1:N
        Zi = cell(B,1);
        for b=1:B,  Zi{b} = blk_ops{b}.Z.in{i};  end
        Zu{i} = unique(vertcat(Zi{:}),'rows','stable');
    end
    [~,~,nZ] = Zdims(Zu);
    % Express each block in the union, with its integrals unique
    for b=1:B
        blk_ops{b} = unique_intmat(blk_ops{b});
        if ~isequal(blk_ops{b}.Z.in,Zu)
            zmap = Zmap(blk_ops{b}.Z.in,Zu);
            P_b = blk_ops{b}.params;
            for k=1:numel(P_b),     P_b{k} = remap_cols(P_b{k},zmap,nZ);    end
            blk_ops{b}.params = P_b;
        end
    end
    % Union of the integrals along each variable: row r of block b is row
    % tgt{b,i}(r) of the union (sorted, unique)
    intmat_new = cell(1,N);     tgt = cell(B,N);        nt = ones(1,max(N,2));
    for i=1:N
        imat_i = cell(B,1);
        for b=1:B,  imat_i{b} = blk_ops{b}.intmat{i};  end
        nr = cellfun(@(a) size(a,1),imat_i);
        [Pmat,intmat_new{i}] = uniquerows_integerTable(vertcat(imat_i{:}));
        [~,tg] = max(Pmat,[],2);
        tg = full(tg);
        off = [0;cumsum(nr)];
        for b=1:B,  tgt{b,i} = tg(off(b)+1:off(b+1));   end
        nt(i) = size(intmat_new{i},1);
    end
    % Add the parameters of each block to its cells in the union. Cells that
    % no block contains keep one shared empty matrix
    Z0 = sparse(size(blk_ops{1}.params{1},1),size(blk_ops{1}.params{1},2));
    params = repmat({Z0},nt);
    is_set = false(nt);
    for b=1:B
        % Linear index in the union of each cell of block b, in the order
        % of blk_ops{b}.params(:) (first variable fastest)
        lin = 1;    stride = 1;
        for i=1:N
            shp = ones(1,max(N,2));     shp(i) = numel(tgt{b,i});
            lin = lin + reshape((tgt{b,i}-1)*stride,shp);
            stride = stride*nt(i);
        end
        for k=1:numel(lin)
            t = lin(k);
            if is_set(t)
                params{t} = params{t} + blk_ops{b}.params{k};
            else
                params{t} = blk_ops{b}.params{k};
                is_set(t) = true;
            end
        end
    end
    Kop = blk_ops{1};
    Kop.params = params;
    Kop.intmat = intmat_new;
    Kop.Z.in = Zu;
    Kop.chkval = true;
end


%% Step 3: Write the kernels of multipliers in terms of the variable they follow
% In a cell, a variable and the multipliers evaluated at it are replaced in
% the kernel by the lowest of them, which is then the one that is
% integrated, and the cell is written accordingly. The functional does not
% change, as apply_op substitutes a multiplier variable by the variable it
% follows. This is done once, after the copies are ordered, as the order
% does not depend on which of the equal variables is integrated.
is_chg = false;
for i=1:N
    imat = Kop.intmat{i};
    nd = size(imat,2)/2;
    if nd<2
        continue
    end
    L = imat(:,1:2:end);        U = imat(:,2:2:end);
    is_dlt = L==U;
    rows_m = find(any(is_dlt,2) & is_cell(imat));
    if isempty(rows_m)
        continue
    end
    % New rows: the variable replacing each variable (the lowest variable
    % of its group) is integrated, the others are evaluated at it
    nm = numel(rows_m);
    nr = zeros(nm,nd);          rows_new = imat(rows_m,:);
    for iq=1:nm
        r = rows_m(iq);
        tgt = 1:nd;             tgt(is_dlt(r,:)) = L(r,is_dlt(r,:));
        rep = accumarray(tgt(:),(1:nd)',[nd,1],@min);
        nr(iq,:) = rep(tgt)';
        for j=1:nd
            if nr(iq,j)==j
                up = U(r,tgt(j));
                if up<=nd,  up = nr(iq,up);     end
                rows_new(iq,2*j-1:2*j) = [0,up];
            else
                rows_new(iq,2*j-1:2*j) = nr(iq,j)*[1,1];
            end
        end
    end
    % Monomials in the variables that are replaced, and the monomials that
    % the parameters in the cells of each new row use
    Zi = Kop.Z.in{i};
    [nmon,st,nZ] = Zdims(Kop.Z.in);
    imgs = cell(nm,1);          chg = false(nmon(i),nm);
    for iq=1:nm
        imgs{iq} = full(Zi*sparse(1:nd,nr(iq,:),1,nd,nd));
        is_new = nr(iq,:)~=(1:nd);
        chg(:,iq) = any(Zi(:,is_new)>0,2);
    end
    used = false(nmon(i),nm);
    sz = size(Kop.params,1:N);
    zi = mod(floor((0:nZ-1)'/st(i)),nmon(i))+1;
    % Row of intmat{i} of each cell (the first variable changes fastest)
    row_cell = mod(floor((0:numel(Kop.params)-1)/prod(sz(1:i-1))),sz(i))+1;
    for iq=find(any(chg,1))
        for lin=find(row_cell==rows_m(iq))
            col_nz = full(any(Kop.params{lin},1));
            used_z = any(reshape(col_nz,nZ,[]),2);
            used_q = any(any(reshape(used_z,st(i),nmon(i),[]),1),3)';
            used(:,iq) = used(:,iq) | (used_q & chg(:,iq));
        end
    end
    need = any(used,1);
    if any(need)
        % Express the monomials in the new basis: the old monomials, and
        % the images of the monomials that are replaced
        Znew = Zi;
        for iq=find(need)
            Znew = [Znew;imgs{iq}(used(:,iq),:)];   %#ok<AGROW>
        end
        Znew = unique(Znew,'rows','stable');
        nmT = nmon;             nmT(i) = size(Znew,1);
        stT = ones(1,N);
        for k=N-1:-1:1
            stT(k) = stT(k+1)*nmT(k+1);
        end
        base = zeros(nZ,1);
        for k=setdiff(1:N,i)
            base = base + mod(floor((0:nZ-1)'/st(k)),nmon(k))*stT(k);
        end
        zmap0 = 1 + base + (zi-1)*stT(i);
        zmaps = cell(1,nm);
        for iq=find(need)
            imap = (1:nmon(i))';
            [~,loc] = ismember(imgs{iq}(used(:,iq),:),Znew,'rows');
            imap(used(:,iq)) = loc;
            zmaps{iq} = 1 + base + (imap(zi)-1)*stT(i);
        end
        for lin=1:numel(Kop.params)
            iq = find(rows_m==row_cell(lin));
            if ~isempty(iq) && need(iq)
                zm = zmaps{iq};
            else
                zm = zmap0;
            end
            Kop.params{lin} = remap_cols(Kop.params{lin},zm,prod(nmT));
        end
        Kop.Z.in{i} = Znew;
        is_chg = true;
    end
    if any(any(rows_new~=imat(rows_m,:)))
        imat(rows_m,:) = rows_new;
        Kop.intmat{i} = imat;
        is_chg = true;
    end
end
if is_chg
    % Add the parameters of rows that now coincide, remove unused monomials
    Kop.chkval = true;
    Kop = unique_intmat(Kop);
    Kop = prune_basis(Kop);
end


%% Step 4: Combine the columns of equal copies of a vector-valued state
% If copies of the same state are equal in all spatial variables
% (multipliers), x_i(s)*x_j(s) = x_j(s)*x_i(s), and only the sum of the
% kernels of these elements is determined: it is placed on the column of
% the product with nondecreasing indices.
has_m = false;
for i=1:N
    imat = Kop.intmat{i};
    has_m = has_m || any(any(imat(:,1:2:end)==imat(:,2:2:end)));
end
if has_m && any(dims_full>1)
    nc = prod(dims_full);
    stc = ones(1,d);
    for j=d-1:-1:1
        stc(j) = stc(j+1)*dims_full(j+1);
    end
    % digs(c,j): element of factor j in column c of the parameters
    digs = mod(floor((0:nc-1)'./stc),dims_full);
    [~,~,nZ] = Zdims(Kop.Z.in);
    sz = size(Kop.params,1:N);      sb = cell(1,N);
    for lin=1:numel(Kop.params)
        [sb{:}] = ind2sub([sz,1],lin);
        % Variable that each factor is equal to, in each spatial variable
        gid = zeros(d,N);
        for i=1:N
            row = Kop.intmat{i}(sb{i},:);
            L = row(1:2:end);       U = row(2:2:end);
            for j=find(depmat_in(:,i))'
                v = pvar_deg(j,i);
                gid(j,i) = v + (L(v)-v)*(L(v)==U(v));
            end
        end
        % Factors of the same state that are equal in all spatial variables
        digs2 = digs;
        for k=unique(state_idcs)
            pos = find(state_idcs==k);
            if numel(pos)<2 || dims_full(pos(1))==1
                continue
            end
            [~,~,grp] = unique(gid(pos,:),'rows');
            for g=1:max(grp)
                pg = pos(grp==g);
                if numel(pg)>1
                    digs2(:,pg) = sort(digs(:,pg),2);
                end
            end
        end
        repcol = 1 + digs2*stc';
        P = Kop.params{lin};
        jcol = find(any(P,1));
        jb = ceil(jcol/nZ);
        if all(repcol(jb)==jb)
            continue
        end
        [ir,jc,vc] = find(P);
        % Column vectors, also if P has one row
        ir = ir(:);     jc = jc(:);     vc = vc(:);
        jb = ceil(jc/nZ);       zc = jc - (jb-1)*nZ;
        jn = (repcol(jb)-1)*nZ+zc;
        Kop.params{lin} = sparse(ir,jn,vc,size(P,1),size(P,2));
    end
end


%% Step 5: Sort the rows (cells) as in unique_intmat
% so that equal functionals have equal intmat.
for i=1:N
    imat = Kop.intmat{i};
    if issortedrows([sum(imat,2),imat])
        continue
    end
    [P,rows_u] = uniquerows_integerTable(imat);
    if size(rows_u,1)==size(imat,1) && ~isequal(rows_u,imat)
        [~,newpos] = max(P,[],2);
        ord = zeros(size(newpos));      ord(newpos) = 1:numel(newpos);
        ip = cellfun(@(a) (1:size(a,1))',Kop.intmat,'UniformOutput',false);
        ip{i} = ord;
        Kop.params = Kop.params(ip{:});
        Kop.intmat{i} = rows_u;
    end
end
end


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function is_ok = is_cell(imat)
% IS_OK = IS_CELL(IMAT) checks which rows of an intmat define an integral
% that is a cell in canonical form: the variables that are not multipliers
% form a chain, each integrated from a up to the next, the largest over
% the full domain, and each multiplier is evaluated at one of these
% variables.
% Used both to find the integrals to expand into cells and the cells whose
% multipliers are rewritten.
d = size(imat,2)/2;
L = imat(:,1:2:end);        U = imat(:,2:2:end);
is_dlt = L==U;
% One variable over [a,b], the others from a; the upper limits are distinct
% variables that are not multipliers; multipliers follow such a variable
is_ok = all(L==0 | is_dlt,2) & sum(U==d+1 & ~is_dlt,2)==1;
cnt = zeros(size(imat,1),d);
for v=1:d
    cnt(:,v) = sum(U==v & ~is_dlt,2);
end
is_ok = is_ok & all(cnt<=~is_dlt,2);
rr = repmat((1:size(imat,1))',1,d);
tgt_dlt = is_dlt(sub2ind(size(is_dlt),rr,max(L,1)));
is_ok = is_ok & all(~is_dlt | ~tgt_dlt,2);
end