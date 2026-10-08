function [A,info] = lift_copvar(dims,spaces,dom,dlift,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [A,INFO] = LIFT_COPVAR(DIMS,SPACES,DOM,DLIFT,OPTIONS) returns the FIXED
% moment lift A of the separated positive-operator form
%
%       Pop = A' * M * A,      M(theta) >= 0 a polynomial multiplier,
%
% as a 'copvar' from the container space X = L_2^{m_1}[s^1] x ... x
% L_2^{m_M}[s^M] (an empty space is R^{m_k}) to the lifted space
%
%       Y = R^{m_R} x L_2^{K_1}[s] x ... x L_2^{K_G}[s],
%
% s the whole registry. The components of Y are the basis operators Z_alpha
% of 'copquadvar' (Sec. 9.1 of the sopvar document) WITHOUT their monomials
% in the auxiliary variable theta: for space k and a multi-index alpha over
% its own variables,
%
%       (A_{k,alpha} x)(theta) = int_{dom(s^k)} I_alpha(theta-s) ...
%                                   (V^alpha(s) kron I_{m_k}) x(s) ds,
%       V^alpha(s) = kron_{d : alpha_d ~= 1} [1; s_d; ...; s_d^{dlift_d}],
%
%       I_1 = delta (multiplier), I_2 lower, I_3 upper, I_4 full integral,
%
% constant in the registry variables s^k lacks, and the identity on R^{m_k}.
% The theta-dependence of copquadvar's basis, Z^alpha(theta,s) =
% Z_w(theta) kron V^alpha(s), moves into M ('posmult_cdopvar'), so that
% A'*M*A with M = g Lambda'*Q*Lambda, Lambda = Z_w(theta) kron I, is the
% quadratic form of 'copquadvar' at degrees int = w, mult = dlift with no
% joint cap. In the proof program's canonical coordinates (proof_roadmap.md
% Sec. 1) the alpha = 2 components are the running moments
% int_a^s V_D(t) x(t) dt and alpha = 4 the full moments, D = dlift the lift
% degree; the alpha = 1 components make it the mixed lift of round 4.
%
% INPUT
% - dims, spaces, dom: as 'copquadvar' (one set of spaces, in = out);
% - dlift:   lift degree per registry variable: a scalar, a 1 x nv row, or
%            a cell with one such entry per space, each of which may be a
%            cell with one entry per component of that space (order of
%            INFO.basis_list). A multiplier direction carries no lift
%            degree. Defaults to 1;
% - options: (optional) struct with fields
%   sep       logical scalar or 1 x nv over the registry: where true the
%             lower and upper components (alpha 2, 3) of that direction are
%             replaced by the full integral (alpha 4), as 'copquadvar';
%   include   1 x M cell, entry k selecting the components of space k in
%             the formats 'copquadvar' accepts: multi-index rows over the
%             space's own variables, a logical mask, or linear indices into
%             the enumeration of 'multiindex_grid'. Defaults to all;
%   ygroup    labels grouping the L2 components into the L2 spaces of Y:
%             a 1 x N_b row (order of INFO.basis_list) or a 1 x M cell with
%             a scalar or a 1 x nb(k) row per space. Components of one group
%             share a weight degree in 'poscopvar_lift'. Labels of R
%             components are ignored. Defaults to one group;
%
% OUTPUT
% - A:       'copvar' with (hasR + G) x M blocks, row 1 the R space of Y
%            when X has an R^{m_k} space (hasR = 1), then one row per
%            group; passes 'verify';
% - info:    struct with fields
%   basis_list  N_b x (1+nv), row c = [space, alpha over the registry, 0
%               where the space lacks the variable]: the component order,
%               that of 'copquadvar' for the same inputs;
%   comp        1 x N_b struct array: space, alpha (own variables), dim,
%               lift (1 x nv), yspace (block row of A), off (row offset in
%               that space, 0-based), ZR (own-variable exponent cells),
%               cells (nonzero parameter cells);
%   dims_Y, spaces_Y, isR_Y  the lifted space as 'poscopvar' takes it, and
%               which of its spaces is R^{m_R};
%   vars, dom, nv, M, sep, nb.
%
% NOTES
% Every component block is canonical by construction: no input degree in a
% multiplier direction, and its only nonzero parameter cells are those the
% indicator selects, gamma_d = alpha_d for alpha_d in {1,2,3} and both
% gamma_d = 2 and 3 for alpha_d = 4. With the output components ordered as
% the ZR monomials (first variable slowest) every nonzero cell is a 0/1
% placement matrix, so no coefficient arithmetic is done here. A block of A
% from space k into a group holds all components of k in that group, each
% on its own rows and on the union of their monomial bases.
%
% Cost: O(N_b) block constructions on the spatial and monomial axes; no
% decision variables. The work of the separated form is the composition
% A'*M*A, whose cost is the block algebra of '@cdopvar/mtimes'.
%
% N-D: the tensor-product structure is the 1-D lift per direction. The
% proof program's results are 1-D, so in N-D this is the extrapolation of
% that structure, not a theorem.
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
    vals = cell(1,nk);
    for t = 1:nk
        if sep(own{k}(t)),  vals{t} = [1,4];  else,  vals{t} = [1,2,3];  end
    end
    Ak = alpha_select(multiindex_grid(vals),incl{k},nk,k);
    if size(Ak,1)==0
        error("At least one component must be included for space "+num2str(k)+".")
    end
    if size(unique(Ak,'rows'),1)~=size(Ak,1)
        error("Multi-indices for space "+num2str(k)+" should be distinct.")
    end
    alpha_sp{k} = Ak;   nb(k) = size(Ak,1);
end
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
        comp(c).dim = dims(k);  comp(c).ZR = cell(1,0);  comp(c).cells = 1;
        comp(c).yspace = 1;
    else
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
        comp(c).ZR = ZR;
        comp(c).dim = dims(k)*prod([cellfun(@numel,ZR),1]);
        comp(c).cells = cell_of_gamma(multiindex_grid(gsets));
        comp(c).yspace = hasR + gl(cL==c);
    end
    comp(c).off = dimY(comp(c).yspace);
    dimY(comp(c).yspace) = dimY(comp(c).yspace) + comp(c).dim;
end

% % % The blocks of A. Row y of Y, column k of X: all components of space k
% in that Y space, on their own rows.
C = cell(nY,M);
for k = 1:M
    if isR(k)
        c = base(k)+1;
        E = sparse(comp(c).off+(1:dims(k)),1:dims(k),1,dimY(1),dims(k));
        C{1,k} = sopvar({E},struct('out',{cell(1,0)},'in',{cell(1,0)}),...
                        cell(1,0),cell(1,0),...
                        struct('out',zeros(0,2),'in',zeros(0,2)),[dimY(1),dims(k)]);
    else
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
Eu = multiindex_grid(ZRu,'first_slowest');      NRu = size(Eu,1);
shape = [3*ones(1,nk),1,1];
params = repmat({sparse(Ky,mk*NRu)},shape);
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
b = sopvar(params,bvars,ZL,ZRu,bdom,[Ky,mk]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = per_component(spec,M,nb,nv,name)
% One 1 x nv row per component from a scalar, a 1 x nv row, a 1 x M cell of
% those, or a 1 x M cell whose entry k is a 1 x nb(k) cell of those.
Nb = sum(nb);   base = cumsum([0,nb]);
out = cell(1,Nb);
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
if isempty(incl)
    return
end
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
