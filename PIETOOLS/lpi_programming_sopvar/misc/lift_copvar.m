function [A,info] = lift_copvar(dims,spaces,dom,dlift,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A,INFO] = LIFT_COPVAR(DIMS,SPACES,DOM,DLIFT,OPTIONS) returns the moment
% lift A, a 'copvar' object without decision variables, from
% X = X_1 x ... x X_M, X_k = R^{m_k} or L_2^{m_k}[s^k], to
%
%       Y = R^{m_R} x L_2^{K_1}[s] x ... x L_2^{K_G}[s],
%
% s the registry, the sorted union of the variables of the s^k. A has one
% component per space k and multi-index alpha over the variables of s^k,
% alpha_d in {1,2,3}, or {1,4} in a separable direction: on L_2^{m_k}[s^k]
%
%       (A_{k,alpha} x)(theta) = int I_alpha(theta-s)
%                                    (I_{m_k} kron V^alpha(s)) x(s) ds,
%       I_alpha(r) = prod_d I_{alpha_d}(r_d),     I_1(r) = delta(r),
%       I_2(r) = [r >= 0],     I_3(r) = [r <= 0],     I_4(r) = 1,
%       V^alpha(s) = kron_{d: alpha_d ~= 1} [1; s_d; ...; s_d^{D_d}],
%
% the integral over the domain of s^k, [.] the indicator function, d over
% the variables of s^k in registry order and D = DLIFT the lift degree;
% the component is constant in the registry variables that s^k lacks. On
% R^{m_k} the one component is the identity. The components of the R^{m_k}
% spaces are stacked in R^{m_R}, and those of the L_2 spaces with one label
% of OPTIONS.ygroup in one L_2^{K_g}[s], each in the order of
% INFO.basis_list; G is the number of labels.
%
% INPUT
% - dims:    'double', M x 1 positive integers, dims(k) = m_k; a scalar
%            for every space;
% - spaces:  1 x M 'cell' of 'cellstr', spaces{k} the variable names of
%            s^k, an empty entry the space R^{m_k}, no name ending in
%            '_int' or '_dum'; a 'cellstr' is one space;
% - dom:     'double', nv x 2, row d the interval of variable d of the
%            registry; or 1 x 2, the interval of every variable; or a
%            'struct' with fields vars ('cellstr') and dom ('double', one
%            row per name);
% - dlift:   (optional) 'double', D: a scalar or a 1 x nv row; or a 1 x M
%            'cell' with one such entry per space, each of which may be a
%            1 x nb(k) 'cell' with one entry per component of the space, in
%            the order of INFO.basis_list; nonnegative integers; default 1;
% - options: (optional) 'struct' with the optional fields
%              sep      'logical', scalar or 1 x nv (default false), the
%                       separable directions;
%              include  1 x M 'cell', entry k the components of space k
%                       (default all, in the enumeration of every alpha,
%                       direction 1 fastest): a 'double' array with one
%                       column per variable of s^k, its rows distinct
%                       admissible alpha; a 'logical' mask over the
%                       enumeration; or other 'double' indices into it;
%              ygroup   'double', the labels of the L_2 components
%                       (default one label): a scalar; or a 1 x M 'cell'
%                       with per space a scalar or a 1 x nb(k) 'cell' of
%                       scalars, one per component; nonnegative integers;
%                       the labels of R^{m_k} components are ignored.
%
% OUTPUT
% - A:     'copvar' object, (hasR + G) x M blocks: block row 1 R^{m_R}
%          when X has an R^{m_k} space (hasR = 1), then one block row per
%          L_2^{K_g}[s]; block (y,k) an 'sopvar' object with the
%          components of space k in space y of Y, [] (zero) where there
%          are none;
% - info:  'struct' with fields
%            basis_list  'double', N_b x (1+nv), row c = [k, alpha] of
%                        component c, alpha over the registry with 0 where
%                        s^k lacks the variable, ordered by space;
%            comp        1 x N_b 'struct' array, per component the fields
%                        space (k), alpha (over the variables of s^k), dim
%                        (its number of rows of Y), lift (1 x nv, D),
%                        yspace (its block row of A), off (its first row
%                        in that space, minus 1), ZR (1 x n_k 'cell', the
%                        exponents of V^alpha in each variable of s^k) and
%                        cells (the parameter cells it fills in its block
%                        of A);
%            dims_Y      'double', (hasR + G) x 1, the dimension of each
%                        space of Y;
%            spaces_Y    1 x (hasR + G) 'cell', the variable names of each
%                        space of Y (an empty 'cell' for R^{m_R});
%            isR_Y       'logical' row, isR_Y(1) true when Y has R^{m_R},
%                        every other entry false;
%            vars, dom   1 x nv 'cellstr', the registry, and nv x 2
%                        'double', its intervals;
%            nv, M, nb   'double', the numbers of registry variables, of
%                        spaces of X and of components per space (1 x M);
%            sep         1 x nv 'logical', the separable directions.
%
% See also POSCOPVAR_LIFT, POSMULT_CDOPVAR, COPQUADVAR, POSCOPVAR,
% MULTIINDEX_GRID, CELL_OF_GAMMA, PARSE_COPVAR_SPACES.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lift_copvar
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
% Initial coding MMP, 10/07/2026. The lift of the separated form A'*M*A
%                proposed from the proof program (its canonical primitive
%                coordinates): the basis operators of 'copquadvar' with the
%                theta-monomials removed, so that the positive multiplier
%                M carries the whole weight. The component enumeration and
%                'include' handling are those of copquadvar (basis_operators,
%                alpha_select), repeated here so that the two routes agree
%                component by component.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was the relation to 'copquadvar' and to the
%                proof program (with a notes reference), implementation
%                notes (canonical blocks, placement matrices) and remarks
%                on cost and N-D. Two help statements now follow the
%                code: ygroup is a scalar or cells of scalars (a row was
%                listed, which the code rejects), and the kernel is
%                I_{m_k} kron V^alpha(s), component outer (V^alpha kron
%                I_{m_k} was written).

% Three inputs at least; DLIFT 1 and OPTIONS struct() by default
if nargin<3
    error("Not enough input arguments.")
end
if nargin<4 || isempty(dlift),      dlift = 1;          end
if nargin<5 || isempty(options),    options = struct(); end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end

% % % Spaces and registry, by the parser 'copquadvar' uses. The registry
% comes back sorted, so every space and every pair's shared set is sorted,
% which is the order a block indexes its parameter cell in.
[meta,spaces,sp_in] = parse_copvar_spaces(dims,spaces,dom);
if ~isequal(sp_in,spaces) || ~isequal(meta.dim_in,meta.dim_out)
    error("A self-adjoint operator has equal input and output spaces.")
end
% Spaces, registry, mask of the variables of each space, component counts,
% domains; the variables of each space as registry indices
M = numel(spaces);      vars = meta.vars;       nv = numel(vars);
mask = meta.space_out;  dims = meta.dim_out(:); dom = meta.dom;
own = cell(1,M);
for k = 1:M
    % Row orientation forced: 'find' on a 1x1 mask returns a column, and a
    % 0x1 empty list breaks the concatenations below.
    own{k} = reshape(find(mask(k,:)),1,[]);
end
isR = cellfun(@isempty,own);

% % % Options
% sep: one flag per registry variable; include: one entry per space
sep = false(1,nv);
if isfield(options,'sep') && ~isempty(options.sep)
    sep = logical(reshape(options.sep,1,[]));
    if isscalar(sep)
        sep = repmat(sep,1,nv);
    elseif numel(sep)~=nv
        error("'sep' should be a scalar or have one entry per registry variable.")
    end
end
incl = cell(1,M);
if isfield(options,'include') && ~isempty(options.include)
    incl = options.include;
    if ~iscell(incl) || iscellstr(incl),    incl = {incl};  end
    incl = reshape(incl,1,[]);
    if numel(incl)~=M
        error("A cell 'include' should have one entry per space.")
    end
end

% % % Components: the multi-indices of copquadvar's basis operators, one
% list per space, direction 1 fastest, then the space's 'include'.
alpha_sp = cell(1,M);   nb = zeros(1,M);
for k = 1:M
    nk = numel(own{k});
    % Admissible entries per variable of the space: 1, 2, 3, or 1, 4 where
    % separable
    vals = cell(1,nk);
    for t = 1:nk
        if sep(own{k}(t)),  vals{t} = [1,4];  else,  vals{t} = [1,2,3];  end
    end
    % The enumeration, then the selection of the space; at least one row,
    % rows distinct
    Ak = alpha_select(multiindex_grid(vals),incl{k},nk,k);
    if size(Ak,1)==0
        error("At least one component must be included for space "+num2str(k)+".")
    end
    if size(unique(Ak,'rows'),1)~=size(Ak,1)
        error("Multi-indices for space "+num2str(k)+" should be distinct.")
    end
    alpha_sp{k} = Ak;   nb(k) = size(Ak,1);
end
% basis_list: row c = [space, alpha over the registry, 0 elsewhere]
Nb = sum(nb);   base = cumsum([0,nb]);
basis_list = zeros(Nb,1+nv);
for k = 1:M
    for i = 1:nb(k)
        c = base(k)+i;
        basis_list(c,1) = k;
        basis_list(c,1+own{k}) = alpha_sp{k}(i,:);
    end
end

% % % Lift degrees, one 1 x nv row per component.
dl = per_component(dlift,M,nb,nv,'dlift');

% % % Groups of Y: the R space (every R component), then one L2 space per
% distinct label among the L2 components, in order of first appearance.
yg = ones(1,Nb);
if isfield(options,'ygroup') && ~isempty(options.ygroup)
    yg = per_component(options.ygroup,M,nb,1,'ygroup');
    yg = cellfun(@(v) v(1),yg);
end
% R and L2 components; G groups of L2 components, numbered in order of
% first appearance
cR = find(isR(basis_list(:,1)));    cR = reshape(cR,1,[]);
cL = find(~isR(basis_list(:,1)));   cL = reshape(cL,1,[]);
hasR = ~isempty(cR);
[~,~,gl] = unique(yg(cL),'stable');     gl = reshape(gl,1,[]);
G = max([gl,0]);
nY = hasR+G;

% % % Per-component geometry: dimension, bases, nonzero cells, Y space and
% row offset.
comp = struct('space',cell(1,Nb),'alpha',[],'dim',[],'lift',[],'yspace',[],...
              'off',[],'ZR',[],'cells',[]);
dimY = zeros(nY,1);
for c = 1:Nb
    k = basis_list(c,1);    ok = own{k};    nk = numel(ok);
    a = alpha_sp{k}(c-base(k),:);
    comp(c).space = k;      comp(c).alpha = a;      comp(c).lift = dl{c};
    if isR(k)
        % R: m_k rows of R^{m_R}, the one parameter cell
        comp(c).dim = dims(k);  comp(c).ZR = cell(1,0);  comp(c).cells = 1;
        comp(c).yspace = 1;
    else
        % L2: per variable of the space, the exponents of V^alpha and the
        % cells of I_alpha
        ZR = cell(1,nk);    gsets = cell(1,nk);
        for t = 1:nk
            switch a(t)
                case 1,     ZR{t} = 0;                  gsets{t} = 1;
                case 2,     ZR{t} = (0:dl{c}(ok(t)))';  gsets{t} = 2;
                case 3,     ZR{t} = (0:dl{c}(ok(t)))';  gsets{t} = 3;
                case 4,     ZR{t} = (0:dl{c}(ok(t)))';  gsets{t} = [2,3];
                otherwise,  error("Internal error: multi-index entry %d.",a(t))
            end
        end
        % Its rows (m_k times the number of monomials of V^alpha), the
        % cells of the product set, the Y space of its group
        comp(c).ZR = ZR;
        comp(c).dim = dims(k)*prod([cellfun(@numel,ZR),1]);
        comp(c).cells = cell_of_gamma(multiindex_grid(gsets));
        comp(c).yspace = hasR + gl(cL==c);
    end
    % Rows off+1..off+dim of its Y space
    comp(c).off = dimY(comp(c).yspace);
    dimY(comp(c).yspace) = dimY(comp(c).yspace) + comp(c).dim;
end

% % % The blocks of A. Row y of Y, column k of X: all components of space k
% in that Y space, on their own rows.
C = cell(nY,M);
for k = 1:M
    if isR(k)
        % R: the identity onto the rows of the component in R^{m_R}
        c = base(k)+1;
        E = sparse(comp(c).off+(1:dims(k)),1:dims(k),1,dimY(1),dims(k));
        C{1,k} = sopvar({E},struct('out',{cell(1,0)},'in',{cell(1,0)}),...
                        cell(1,0),cell(1,0),...
                        struct('out',zeros(0,2),'in',zeros(0,2)),[dimY(1),dims(k)]);
    else
        % L2: one block per Y space that holds components of space k
        for y = (hasR+1):nY
            cs = find(basis_list(:,1)==k & [comp.yspace]'==y);
            if isempty(cs),     continue,   end
            C{y,k} = lift_block(comp(cs),dimY(y),own{k},vars,dom,dims(k));
        end
    end
end

% % % The container, with its metadata given (nothing to re-derive).
so = false(nY,nv);      so((hasR+1):nY,:) = true;
mt = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',mask,...
            'dim_out',dimY,'dim_in',dims);
A = copvar(C,mt);

% The lifted space and the component data
spaces_Y = cell(1,nY);
if hasR,    spaces_Y{1} = cell(1,0);    end
for y = (hasR+1):nY,    spaces_Y{y} = vars;    end
info = struct('basis_list',basis_list,'comp',comp,'dims_Y',dimY,...
              'spaces_Y',{spaces_Y},'isR_Y',[hasR,false(1,G)],...
              'vars',{vars},'dom',dom,'nv',nv,'M',M,'sep',sep,'nb',nb);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function b = lift_block(comps,Ky,ok,vars,dom,mk)
% The block of A from space k (own variables ok, mk components) into a Y
% space of dimension Ky, holding the components COMPS: rows comp.off+1..,
% columns on the union of their monomial bases, each nonzero cell a 0/1
% placement. Output variables [S2,S3] with S2 the registry variables k
% lacks (constant there), input variables S3 = s^k.
% INPUT  comps: 'struct' array, the entries of INFO.comp in the block;
%        Ky:    'double', the dimension of the Y space;
%        ok:    1 x n_k 'double', the registry indices of the variables of
%               space k;
%        vars:  1 x nv 'cellstr', the registry; dom: nv x 2 'double';
%        mk:    'double', the number of components m_k of space k.
% OUTPUT b:     'sopvar' object, Ky x mk.
% Output variables [S2,S3], input variables S3, their domains; a constant
% left basis
nv = numel(vars);   nk = numel(ok);
S3 = vars(ok);      S2 = vars(setdiff(1:nv,ok));
vout = [S2,S3];     vin = S3;
[~,io] = ismember(vout,vars);   [~,ii] = ismember(vin,vars);
bvars = struct('out',{vout},'in',{vin});
bdom  = struct('out',dom(io,:),'in',dom(ii,:));
ZL = repmat({0},1,numel(vout));

% Union right basis, per own direction; a direction that is a multiplier in
% every component stays at degree 0 (canonical multiplier form).
ZRu = cell(1,nk);
for t = 1:nk
    e = 0;
    for c = 1:numel(comps),     e = max(e,max(comps(c).ZR{t}));    end
    ZRu{t} = (0:e)';
end
% Union monomials, first variable slowest; a zero matrix per parameter cell
Eu = multiindex_grid(ZRu,'first_slowest');      NRu = size(Eu,1);
shape = [3*ones(1,nk),1,1];
params = repmat({sparse(Ky,mk*NRu)},shape);
% Each component: its placement matrix E1, added to each of its cells
for c = 1:numel(comps)
    Ec = multiindex_grid(comps(c).ZR,'first_slowest');  Nc = size(Ec,1);
    [tf,pos] = ismember(Ec,Eu,'rows');
    if ~all(tf),    error("Internal error: a component monomial is not in the union basis.");   end
    % Output rows (component outer, monomial inner) against columns
    % (component outer, union monomial inner): entry 1 where the monomials
    % agree, as the kernel is V^alpha(s') kron I_mk.
    [J,I] = meshgrid(1:Nc,1:mk);
    rows = comps(c).off + (I(:)-1)*Nc + J(:);
    cols = (I(:)-1)*NRu + pos(J(:));
    E1 = sparse(rows,cols,1,Ky,mk*NRu);
    for q = reshape(comps(c).cells,1,[])
        params{q} = params{q} + E1;
    end
end
% The block
b = sopvar(params,bvars,ZL,ZRu,bdom,[Ky,mk]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = per_component(spec,M,nb,nv,name)
% One 1 x nv row per component from a scalar, a 1 x nv row, a 1 x M cell of
% those, or a 1 x M cell whose entry k is a 1 x nb(k) cell of those.
% INPUT  spec: the DLIFT or ygroup input;
%        M:    'double', the number of spaces;
%        nb:   1 x M 'double', the number of components of each space;
%        nv:   'double', the length of a row;
%        name: 'char', the name used in the error messages.
% OUTPUT out:  1 x sum(nb) 'cell' of 1 x max(nv,1) 'double' rows.
Nb = sum(nb);   base = cumsum([0,nb]);
out = cell(1,Nb);
% A cell: per space one value, or a cell with one per component;
% otherwise one value for every component
if iscell(spec)
    if numel(spec)~=M
        error("A cell '%s' should have one entry per space.",name)
    end
    for k = 1:M
        sk = spec{k};
        if iscell(sk)
            if numel(sk)~=nb(k)
                error("A cell '%s' for space %d should have one entry per component (%d).",...
                      name,k,nb(k))
            end
            for i = 1:nb(k),    out{base(k)+i} = expand_row(sk{i},nv,name);    end
        else
            for i = 1:nb(k),    out{base(k)+i} = expand_row(sk,nv,name);       end
        end
    end
else
    for c = 1:Nb,   out{c} = expand_row(spec,nv,name);  end
end
end


function v = expand_row(v,nv,name)
% V = EXPAND_ROW(V,NV,NAME) returns V as a row of nonnegative integers,
% one entry per registry variable, a scalar repeated.
% INPUT  v:    'double', a scalar or nv entries;
%        nv:   'double', the number of registry variables;
%        name: 'char', the name used in the error messages.
% OUTPUT v:    1 x max(nv,1) 'double'.
if ~isnumeric(v) || any(v(:)<0) || any(v(:)~=round(v(:)))
    error("'%s' should be nonnegative integers.",name)
end
v = reshape(v,1,[]);
if isscalar(v)
    v = repmat(v,1,max(nv,1));
elseif numel(v)~=nv
    error("'%s' should be a scalar or have one entry per registry variable.",name)
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = alpha_select(A,incl,nd,k)
% One space's 'include', in the three formats 'copquadvar' accepts (its
% local function of the same name): multi-index rows validated against the
% generated list, a logical mask, or linear indices.
% INPUT  A:    'double', the enumeration of the multi-indices, one per row;
%        incl: the entry of OPTIONS.include for the space ([] for all);
%        nd:   'double', the number of variables of the space;
%        k:    'double', the space index, for the error messages.
% OUTPUT A:    'double', the selected rows.
% Empty: every row
if isempty(incl)
    return
end
% A mask, rows validated against the enumeration, or linear indices
if islogical(incl)
    if numel(incl)~=size(A,1)
        error("A logical 'include' for space "+num2str(k)+" has one entry "...
              +"per component of that space.")
    end
    A = A(incl(:),:);
elseif nd>0 && size(incl,2)==nd
    if ~all(ismember(incl,A,'rows'))
        error("'include' for space "+num2str(k)+" contains an inadmissible "...
              +"multi-index.")
    end
    A = incl;
else
    if any(incl(:)<1) || any(incl(:)>size(A,1)) || any(incl(:)~=round(incl(:)))
        error("Linear indices in 'include' for space "+num2str(k)...
              +" are out of range.")
    end
    A = A(incl(:),:);
end
end
