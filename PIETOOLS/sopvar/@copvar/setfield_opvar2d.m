function P = setfield_opvar2d(P,fname,val,cidx,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = setfield_opvar2d(P,fname,val) sets the field fname of the PI operator
%   P: R^m0 x L2^mx[x] x L2^my[y] x L2^m2[x,y]
%                           -> R^n0 x L2^nx[x] x L2^ny[y] x L2^n2[x,y]
% as if P were an 'opvar2d' object, with the components Rij mapping the
% space j to the space i (i,j = 0, x, y, 2), so that P.(fname) then
% returns val.
%
% INPUT
% P:      copvar class object with the structure of an 'opvar2d': at most
%         two spatial variables, and the block rows (columns) ordered by
%         space, first R, then L2[x], then L2[y], then L2[x,y]
% fname:  char or string, the name of an 'opvar2d' field: 'R00', ...,
%         'R22', 'I', 'var1' or 'var2'
% val:    the value of the field:
%         'Rij':    ni x mj 'double' or 'polynomial', or a cell of these,
%                   3 x 1 if x is in both spaces i and j, 1 x 3 if y is,
%                   3 x 3 if both are, [] for zero. Kernels are in var1 (x,y) and, in
%                   a direction with an integral over part of the domain, in
%                   var2 (x_dum, y_dum); a variable only in the input space
%                   may be named x or x_dum;
%         'I':      nv x 2 array, the domain of each of the nv variables;
%         'var1':   nv x 1 'polynomial', the spatial variables [x;y], nv <= 2;
%         'var2':   nv x 1 'polynomial', the dummy variables; only the
%                   names of var1 with '_dum' appended are accepted;
% cidx:   (optional) linear index, or cell of subscripts, of one parameter
%         of a cell component: only that parameter is set, to the kernel
%         val. E.g. {2,3} for P.R22{2,3} = val;
% opts:   (optional) struct with fields
%         vars:     cellstr, names for variables P does not have yet when
%                   the kernel names none, default {'s1','s2'};
%         label:    char, the field name in error messages, default fname;
%
% OUTPUT
% P:      the copvar with the field set
%
% NOTES:
% - x is P.vars{1} and y is P.vars{2} (sorted).
%   A variable P does not have yet is taken from the kernel (its name with
%   '_dum' removed), else from opts.vars, with domain [0,1]; a second
%   variable must sort after the first, else it would be read as x.
% - Only the blocks of the component change, so the block rows of P are
%   kept. A component on a space of dimension 0 adds a block row (column)
%   for it, in the order R, L2[x], L2[y], L2[x,y], the other components in
%   it being zero. A row or column of the component left without a block
%   keeps a zero block, as the container reads its space off a block.
%   Setting a component with no rows or columns does nothing.
% - var1 renames the variables only if no block depends on them, or adds
%   y. I changes the domain of every block.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - setfield_opvar2d (copvar)
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/15/2026: Initial coding;

% Defaults of the optional inputs; check the structure of P and the field
if nargin<4
    cidx = [];
end
if nargin<5
    opts = struct();
end
if ~isfield(opts,'vars')
    opts.vars = {'s1','s2'};
end
if ~isfield(opts,'label')
    opts.label = fname;
end
% P must have the spaces R, L2[x], L2[y], L2[x,y] of an 'opvar2d'
[is_2d,idx_out,idx_in] = is_opvar2d(P);
if ~is_2d
    error('copvar:notOpvar2d',['The operator does not have the 2-D field layout (R00 to R22) of an operator '...
        'R x L2[x] x L2[y] x L2[x,y] -> R x L2[x] x L2[y] x L2[x,y]; see ''is_opvar2d''.'])
end
if ~(ischar(fname) || (isstring(fname) && isscalar(fname)))
    error('copvar:badField','The field name must be a char or string.')
end
fname = char(fname);    label = char(opts.label);
% Component (i,j) maps space j to space i, the spaces 1 to 4 being R,
% L2[x], L2[y], L2[x,y]; has_x(i) (has_y(i)) is true if space i has x (y).
Rnames = {'R00','R0x','R0y','R02';
          'Rx0','Rxx','Rxy','Rx2';
          'Ry0','Ryx','Ryy','Ry2';
          'R20','R2x','R2y','R22'};
has_x = [false,true,false,true];        has_y = [false,false,true,true];
vnames = P.vars;                        % x = vnames{1}, y = vnames{2}

% % The variables and domain, which no component depends on
switch fname
    case 'dim'
        error('copvar:setDim','''dim'' follows from the components; set those instead.')
    case {'var1','var2'}
        % The names of a vector of pvars
        if ~isa(val,'polynomial'),  val = polynomial(val);  end
        vnames_new = cell(1,numel(val));
        for k=1:numel(val)
            v = val(k);
            if numel(v.varname)~=1 || ~isequal(full(v.coefficient),1) || ~isequal(full(v.degmat),1)
                error('copvar:badValue','''%s'' must be a vector of spatial variables (pvar).',label)
            end
            vnames_new(k) = v.varname;
        end
        if strcmp(fname,'var2')
            % A 'sopvar' records no dummy name: only the default is possible
            if isempty(vnames),     vnames = opts.vars(1:min(numel(vnames_new),end));   end
            if ~isequal(vnames_new,strcat(vnames,'_dum'))
                error('copvar:dummyName',['Dummy variables must have the default names, the '...
                    'spatial variables with ''_dum'' appended: %s. Set var1 first to use other names.'],...
                    strjoin(strcat(vnames,'_dum'),', '))
            end
            return
        end
        % var1: other names only if no block depends on the variables, or
        % a variable y after x
        if isequal(vnames_new,vnames),  return,     end
        has_L2 = any(P.space_out(:)) || any(P.space_in(:));
        if numel(vnames_new)>2 || ~issorted(vnames_new)
            error('copvar:varOrder',['var1 must hold at most two variables, in alphabetical order '...
                '(the container stores them sorted and reads the first as x).'])
        elseif has_L2 && ~(isscalar(vnames) && numel(vnames_new)==2 && strcmp(vnames_new{1},vnames{1}))
            error('copvar:renameVars',['The variables of an operator with L2 components cannot be '...
                'renamed; it has {%s}, and var1 may only add a second variable after these.'],strjoin(vnames,','))
        end
        nv = numel(vnames_new);
        P = set_vars(P,vnames_new,[P.dom(1:min(end,nv),:); repmat([0,1],nv-min(size(P.dom,1),nv),1)]);
        return
    case 'I'
        if ~isnumeric(val) || size(val,2)~=2 || any(val(:,1)>=val(:,2))
            error('copvar:badValue','''%s'' must be an nv x 2 array [a,b] of domains, with a < b.',label)
        elseif isempty(vnames) && size(val,1)<=numel(opts.vars)
            % Name the variables of the domain by default
            P = set_vars(P,opts.vars(1:size(val,1)),double(val));
            return
        elseif size(val,1)~=numel(vnames)
            error('copvar:badValue','''%s'' must have one row per spatial variable, %d.',label,numel(vnames))
        end
        % Every block on the new domain; the kernels do not depend on it
        C = P.C;
        for k=find(~cellfun(@isempty,C(:)))'
            blk = C{k};
            [~,idcs_out] = ismember(blk.vars.out,vnames);   [~,idcs_in] = ismember(blk.vars.in,vnames);
            dom = struct('out',val(idcs_out,:),'in',val(idcs_in,:));
            C{k} = sopvar(blk.params,blk.vars,blk.ZL,blk.ZR,dom,blk.dims);
        end
        meta = metadata(P);     meta.dom = double(val);
        P = copvar(C,meta);
        return
end

% % A component Rij; its parameters form a cell of size sz
[i,j] = find(strcmp(Rnames,fname));
if isempty(i)
    error('copvar:badField',['''%s'' is not in the 2-D field layout; use one of '...
        'R00, R0x, R0y, R02, Rx0, Rxx, Rxy, Rx2, Ry0, Ryx, Ryy, Ry2, R20, '...
        'R2x, R2y, R22, I, var1, var2.'],fname)
end
sz = [1+2*(has_x(i) && has_x(j)), 1+2*(has_y(i) && has_y(j))];
% The kernels to set: all parameters, or the one at cidx
if isempty(cidx)
    idcs_K = 1:prod(sz);
    if prod(sz)==1,     K = {val};
    elseif iscell(val) && isequal(size(val),sz),    K = val;
    else
        error('copvar:badValue','''%s'' must be a %d x %d cell of kernels.',label,sz(1),sz(2))
    end
else
    if iscell(cidx) && isscalar(cidx)
        idcs_K = cidx{1};                   % a linear index
    elseif iscell(cidx)
        try
            idcs_K = sub2ind(sz,cidx{:});
        catch
            idcs_K = [];                    % out of range: refused below
        end
    else
        idcs_K = cidx;
    end
    if ~isscalar(idcs_K) || idcs_K<1 || idcs_K>prod(sz)
        error('copvar:badIndex','The index of the parameter of ''%s'' must select one of its %d parameters.',label,prod(sz))
    end
    K = cell(sz);   K{idcs_K} = val;
end

% The size of the component: that of P, or else that of the kernels; []
% is zero of that size
for k=idcs_K
    if isa(K{k},'dpvar')
        error('copvar:decisionValue',['A copvar holds no decision variables; convert it with '...
            'cdopvar(P) to set ''%s'' to a dpvar.'],label)
    elseif ~(isnumeric(K{k}) || islogical(K{k}) || isa(K{k},'polynomial'))
        error('copvar:badValue','''%s'' must be a ''double'' or ''polynomial'' matrix.',label)
    elseif ~isempty(K{k})
        K{k} = polynomial(K{k});
    end
end
% Rows and columns of the component in P, and the one size of the kernels
nr_P = sum(P.dim_out(idx_out==i));      nc_P = sum(P.dim_in(idx_in==j));
sz_K = cellfun(@size,K(idcs_K),'UniformOutput',false);
sz_K = unique(cell2mat(sz_K(~cellfun(@isempty,K(idcs_K)))'),'rows');
if size(sz_K,1)>1
    error('copvar:sizeMismatch','The parameters of ''%s'' must all have the same size.',label)
elseif isempty(sz_K)
    sz_K = [nr_P,nc_P];
end
nr = sz_K(1);   nc = sz_K(2);
% The kernels must match the nonzero dimensions of the spaces of P; a
% component of size zero leaves nothing to set
spaces = {'R','L2[x]','L2[y]','L2[x,y]'};
if nr_P>0 && nr~=nr_P
    error('copvar:sizeMismatch','''%s'' must have %d rows, the dimension of the %s output space; it has %d.',label,nr_P,spaces{i},nr)
elseif nc_P>0 && nc~=nc_P
    error('copvar:sizeMismatch','''%s'' must have %d columns, the dimension of the %s input space; it has %d.',label,nc_P,spaces{j},nc)
elseif nr==0 || nc==0
    return
end
for k=idcs_K
    if isempty(K{k}),   K{k} = polynomial(zeros(nr,nc));    end
end

% The variables x (and y) of the component: those of P, else those the
% kernels name ('_dum' removed), else the defaults
nv = max(1*(has_x(i) || has_x(j)), 2*(has_y(i) || has_y(j)));
if numel(vnames)<nv
    vnames_K = cell(0,1);
    for k=idcs_K
        vnames_K = [vnames_K; K{k}.varname(any(K{k}.degmat,1))]; %#ok<AGROW>
    end
    vnames_new = setdiff(regexprep(vnames_K,'_dum$',''),vnames);
    if isempty(vnames_new) && numel(opts.vars)>=nv && ~any(ismember(opts.vars(numel(vnames)+1:nv),vnames))
        vnames_new = opts.vars(numel(vnames)+1:nv);     % unless a default is taken
    end
    vnames_new = [vnames, reshape(vnames_new,1,[])];
    if numel(vnames_new)~=nv
        error('copvar:needVars',['''%s'' needs %d spatial variable(s), and the operator has {%s}: '...
            'set var1 first, as [x;y] in alphabetical order.'],label,nv,strjoin(vnames,','))
    elseif ~issorted(vnames_new)
        error('copvar:varOrder',['The variable ''%s'' would sort before ''%s'' and be read as x; '...
            'name the variables in alphabetical order, x before y.'],vnames_new{end},vnames_new{1})
    end
    P = set_vars(P,vnames_new,[P.dom; repmat([0,1],nv-size(P.dom,1),1)]);
    vnames = vnames_new;
end
nv = numel(vnames);
vars.out = vnames([has_x(i),has_y(i)] & (1:2)<=nv);
vars.in = vnames([has_x(j),has_y(j)] & (1:2)<=nv);
dom.out = P.dom(ismember(vnames,vars.out),:);   dom.in = P.dom(ismember(vnames,vars.in),:);
% An input variable also in the output is integrated over part of its
% domain, and named with '_dum' in the kernels (vars_K)
is_shared = ismember(vars.in,vars.out);
vars_K = vars;      vars_K.in(is_shared) = strcat(vars.in(is_shared),'_dum');

% Check the variables of each kernel. An input variable integrated over all
% of its domain may be named with '_dum'; a multiplier direction (parameter
% index 1 along it) cannot depend on its dummy variable.
for k=idcs_K
    Kk = K{k};      vnames_k = Kk.varname(any(Kk.degmat,1));
    for v=vars.in(~is_shared)
        is_dum = strcmp(Kk.varname,[v{1},'_dum']);
        if any(is_dum) && ismember(v{1},vnames_k)
            error('copvar:badKernelVar','''%s'' integrates over all of %s, so it cannot depend on both %s and %s_dum.',label,v{1},v{1},v{1})
        elseif any(is_dum)
            vnames_k = Kk.varname;      vnames_k(is_dum) = v;
            Kk = polynomial(Kk.coefficient,Kk.degmat,vnames_k,Kk.matdim);
            vnames_k = Kk.varname(any(Kk.degmat,1));
        end
    end
    vnames_bad = setdiff(vnames_k,[vars_K.out,vars_K.in]);
    if ~isempty(vnames_bad)
        error('copvar:badKernelVar','''%s'' depends on ''%s''; it may only depend on {%s}.',label,vnames_bad{1},strjoin([vars_K.out,vars_K.in],', '))
    end
    gamma_k = [1,1];    [gamma_k(1),gamma_k(2)] = ind2sub(sz,k);    % parameter index along x and y
    dum_mult = vars_K.in(is_shared);
    dum_mult = intersect(dum_mult(gamma_k(sz==3)==1),vnames_k);
    if ~isempty(dum_mult)
        where = sprintf('''%s''',label);
        if isempty(cidx),   where = sprintf('Parameter %d of %s',k,where);    end
        error('copvar:multiplierDummy','%s is a multiplier along %s, so it cannot depend on the dummy variable %s.',...
            where,dum_mult{1}(1:end-4),dum_mult{1})
    end
    K{k} = Kk;
end

% A space of dimension 0 gets a block row (column) of the size of the
% kernels, placed so that the order R, L2[x], L2[y], L2[x,y] is kept
C = P.C;    dim_out = P.dim_out;    dim_in = P.dim_in;
if nr_P==0
    r = sum(idx_out<=i);
    C = [C(1:r,:); cell(1,size(C,2)); C(r+1:end,:)];
    idx_out = [idx_out(1:r); i; idx_out(r+1:end)];      dim_out = [dim_out(1:r); nr; dim_out(r+1:end)];
end
if nc_P==0
    c = sum(idx_in<=j);
    C = [C(:,1:c), cell(size(C,1),1), C(:,c+1:end)];
    idx_in = [idx_in(1:c); j; idx_in(c+1:end)];         dim_in = [dim_in(1:c); nc; dim_in(c+1:end)];
end

% The blocks of the component, each from its part of the kernels
rows = find(idx_out==i);    cols = find(idx_in==j);
sz_params = [3*ones(1,nnz(is_shared)),1,1];         % parameters of a block
K_zero = @(m,n) repmat({polynomial(zeros(m,n))},sz_params);
for a=1:numel(rows)
    for b=1:numel(cols)
        r = rows(a);    c = cols(b);
        ridcs = sum(dim_out(rows(1:a-1))) + (1:dim_out(r));
        cidcs = sum(dim_in(cols(1:b-1))) + (1:dim_in(c));
        K_blk = K_zero(dim_out(r),dim_in(c));
        for k=idcs_K
            K_blk{k} = K{k}(ridcs,cidcs);
        end
        is_zero = all(cellfun(@(Kk) nnz(Kk.coefficient)==0,K_blk(:)));
        blk = C{r,c};
        if ~isempty(blk) && numel(idcs_K)<prod(sz)
            % Keep the other parameters: zero these, add the new ones
            nZL = prod(cellfun(@numel,blk.ZL));     nZR = prod(cellfun(@numel,blk.ZR));
            params = blk.params;
            for k=idcs_K
                params{k} = sparse(dim_out(r)*nZL,dim_in(c)*nZR);
            end
            blk = sopvar(params,blk.vars,blk.ZL,blk.ZR,blk.dom,blk.dims);
            if ~is_zero
                blk = blk + build_block(K_blk,vars,vars_K,dom,[dim_out(r),dim_in(c)]);
            end
        elseif ~is_zero
            blk = build_block(K_blk,vars,vars_K,dom,[dim_out(r),dim_in(c)]);
        else
            blk = [];
        end
        C{r,c} = blk;
    end
end
% A row or column of the component left without a block keeps a zero one
for r=rows'
    if all(cellfun(@isempty,C(r,:)))
        C{r,cols(1)} = build_block(K_zero(dim_out(r),dim_in(cols(1))),vars,vars_K,dom,[dim_out(r),dim_in(cols(1))]);
    end
end
for c=cols'
    if all(cellfun(@isempty,C(:,c)))
        C{rows(1),c} = build_block(K_zero(dim_out(rows(1)),dim_in(c)),vars,vars_K,dom,[dim_out(rows(1)),dim_in(c)]);
    end
end

% The metadata. The one-argument constructor is not called to validate the
% blocks: it refuses a row without blocks, which PIE operators may have.
space_out = [has_x(idx_out)', has_y(idx_out)'];     space_in = [has_x(idx_in)', has_y(idx_in)'];
meta = struct('vars',{vnames},'dom',P.dom,'space_out',space_out(:,1:nv),'space_in',space_in(:,1:nv),...
    'dim_out',dim_out,'dim_in',dim_in);
P = copvar(C,meta);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function blk = build_block(K,vars,vars_K,dom,dims)
% blk = build_block(K,vars,vars_K,dom,dims) is the 'sopvar' block with the
% kernels K.
%
% INPUT
% K:        cell of 'polynomial' kernels, the shape of the block's
%           parameters, each dims(1) x dims(2)
% vars:     struct with fields 'out' and 'in', cellstr of the output and
%           input variables of the block
% vars_K:   struct as vars, with the names of the input variables in the
%           kernels (v_dum for a variable also in the output)
% dom:      struct with fields 'out' and 'in', the domains of vars
% dims:     1 x 2, the dimensions of the block
%
% OUTPUT
% blk:      'sopvar' with blk.kernels equal to K
%
% NOTES
% 'dpvars2sdvars' writes all kernels over common monomial bases, the
% coefficients as vec(C) (component outer, monomial inner).
[A,~,Z] = dpvars2sdvars(K,vars_K);
nZL = prod(cellfun(@numel,Z.out));      nZR = prod(cellfun(@numel,Z.in));
params = cellfun(@(Ak) reshape(Ak,dims(1)*nZL,dims(2)*nZR),A,'UniformOutput',false);
blk = sopvar(params,vars,Z.out,Z.in,dom,dims);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function P = set_vars(P,vnames,dom)
% P = set_vars(P,vnames,dom) gives the copvar P the spatial variables
% vnames, on the domains dom, keeping its blocks.
%
% INPUT
% P:        copvar, whose blocks depend on at most the first variables of
%           vnames, as in P.vars
% vnames:   1 x nv cellstr, sorted
% dom:      nv x 2 array
%
% OUTPUT
% P:        the copvar with P.vars = vnames and P.dom = dom
%
% NOTES
% The masks of the spaces get a column of false per new variable;
% A(:,end+1) = false would instead add a row to an empty (0 x 0) mask.
meta = metadata(P);     nv = numel(vnames);
meta.vars = vnames;     meta.dom = dom;
meta.space_out = [meta.space_out(:,1:min(end,nv)), false(numel(P.dim_out),nv-min(size(meta.space_out,2),nv))];
meta.space_in = [meta.space_in(:,1:min(end,nv)), false(numel(P.dim_in),nv-min(size(meta.space_in,2),nv))];
P = copvar(P.C,meta);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
