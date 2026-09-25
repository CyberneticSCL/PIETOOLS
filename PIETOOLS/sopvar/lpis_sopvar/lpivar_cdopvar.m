function [prog,Pop] = lpivar_cdopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP] = LPIVAR_CDOPVAR(PROG,DIMS,SPACES,DOM,DEG,OPTIONS) declares a
% GENERAL decision operator between concatenations of mixed spaces,
%
%   Pop: L_2^{p_1}[s^1] x ... x L_2^{p_N}[s^N]
%                   -> L_2^{q_1}[t^1] x ... x L_2^{q_M}[t^M],
%
% whose every kernel coefficient is a free decision variable: no positivity,
% no self-adjointness, and the input and output spaces may differ. It is the
% container counterpart of 'lpivar' (and 'lpivar_2d'), which the Q-form
% stability and H-infinity executives use for Qop, and which control and
% estimation use for the rectangular gain variable Zop.
%
% Each block (i,j) is an 'sdopvar' from space j to space i. With S3 the
% variables the two spaces share, S2 those only in the output and S1 those
% only in the input, its kernel in gamma cell gam (Sec. 4-5) is
%
%   K_gam(t,s') = (I_q (x) ZL(t))' C_gam (I_p (x) ZR(s')),
%
% gam_k = 1 a multiplier delta(t_k - s'_k), 2 and 3 the lower and upper
% integrals in direction k of S3. Every admissible entry of every C_gam is
% its own decision variable; admissible means within the degree caps below
% and, in a multiplier direction, at right degree 0, which is the canonical
% multiplier form the class requires (Sec. 5), so the constructor has
% nothing to fold.
%
% INPUT
% - prog:   'struct' specifying an LPI program to modify;
% - dims:   struct with fields 'out' (M x 1) and 'in' (N x 1), the component
%           counts q_i and p_j; or an M x 1 array for a square container
%           whose input spaces are its output spaces. A scalar is expanded;
% - spaces: struct with fields 'out' (1 x M cell) and 'in' (1 x N cell) of
%           'cellstr' variable names, one per space, an empty entry being
%           the finite-dimensional space R^q; or a 1 x M cell for a square
%           container, as 'poscopvar' takes it. A plain 'cellstr' is one
%           space;
% - dom:    nv x 2 domains in the order of the SORTED registry of all
%           variables, a single 1 x 2 row for every variable, or a struct
%           with fields 'vars' and 'dom' pairing names with intervals - the
%           forms 'poscopvar' accepts;
% - deg:    (optional, default 1) the degree caps, per variable ROLE:
%             'mult'  left degree in a multiplier direction of S3;
%             'int'   [left, right] degrees in an integral direction of S3;
%             'out'   degree in an output-only (S2) variable;
%             'in'    degree in an input-only (S1) variable.
%           A scalar d sets mult = out = in = d and int = [d d]. A 1 x 3
%           array [d1 d2 d3] is 'lpivar's convention and sets mult = out =
%           in = d1, int = [d3 d2]: in 1-D that reproduces lpivar's P, Q1,
%           Q2, R0 (degree d1) and R1, R2 (d3 in s, d2 in theta). A struct
%           may give any subset of the fields, the rest defaulting to 1. An
%           M x N cell gives one specification per block;
% - options: (optional) struct with field
%           - occ:  M x N logical, false for a structurally zero block,
%                   which declares no variables. Default all true;
%                   a row or column entirely false still holds one          % MMP, 09/25/2026
%                   explicit zero block, with no variables;                 % MMP, 09/25/2026
%
% OUTPUT
% - prog:   the program with the new decision variables declared;
% - Pop:    'cdopvar' object. All blocks are on one decision variable list,
%           declared by a single 'lpidecvar' call.
%
% NOTES
% Every decision variable multiplies exactly one monomial entry of one
% kernel, so the number of variables equals the dimension of the operator
% family and the map from variables to operators is injective. For 1-D
% spaces and lpivar's degrees the count is lpivar's own,
%   n11*n12 + (d1+1)*(n11*n22 + n21*n12 + n21*n22) + 2*(d3+1)*(d2+1)*n21*n22,
% and the kernels span the same family (see 'test_lpivar_cdopvar').
%
% Degree caps are per ROLE, not per variable. Unlike 'lpivar_2d' there is
% no 'sep' option: an operator separable in a direction is a restriction of
% this family, and can be imposed with 'lpi_eq_cdopvar' if wanted.
%
% The domain parsing mirrors 'copquadvar', which does it inline, rather than
% refactor that routine into a shared helper.
%
% Cost: O(number of decision variables) - one triplet per variable, one
% 'sparse' per gamma cell, and one declaration for the whole container. The
% decision rows of every block's B are the container's full list, which
% costs only B's column pointers, O(nC), not O(q) per block.
%
% See also LPIVAR, LPIVAR_2D, POSCOPVAR, LPI_EQ_CDOPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpivar_cdopvar
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
% Initial coding MMP, 09/25/2026
% MMP, 09/25/2026: A row or column that 'occ' leaves with no block gets one
%                  explicit zero block (degree 0, no variables). Without it
%                  the container failed 'verify' and could not be rebuilt by
%                  cdopvar(C), which read a row's space off a block. Costs
%                  O(nC) of column pointers for that block, nothing in q.

if nargin<5 || isempty(deg),        deg = 1;            end
if nargin<6 || isempty(options),    options = struct(); end
if ~isa(prog,'struct') || ~isfield(prog,'decvartable')
    error("'prog' should be an LPI program structure; use 'lpiprogram'.")
end

% % % Spaces and dimensions
[sp_out,sp_in,d_out,d_in] = parse_spaces(dims,spaces);
M = numel(sp_out);      N = numel(sp_in);

% % % Registry, sorted, as the containers and 'poscopvar' keep it: an
% sdopvar indexes its gamma cells over S3 in sorted order, and with a sorted
% registry every space and shared set is already sorted.
vars = reshape(unique([sp_out{:}, sp_in{:}]),1,[]);     nv = numel(vars);
is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));
if any(is_reserved)
    error("Spatial variable names may not end in '_int' or '_dum'; "...
          +"'"+string(vars{find(is_reserved,1)})+"' does.")
end
dom = parse_dom(dom,vars);

% % % Occupancy and degrees
occ = true(M,N);
if isfield(options,'occ') && ~isempty(options.occ)
    occ = logical(options.occ);
    if ~isequal(size(occ),[M,N])
        error("options.occ should be %d x %d, one entry per block.",M,N)
    end
end
degs = parse_deg(deg,M,N);
% A row or column with no block cannot have its space read off one          % MMP, 09/25/2026
% ('verify', cdopvar(C)), so it gets one explicit zero block: degree 0,     % MMP, 09/25/2026
% no coefficient positions, hence no variables. Rows first, then columns.   % MMP, 09/25/2026
zfill = false(M,N);                                                         % MMP, 09/25/2026
for i = find(~any(occ,2))',         zfill(i,1) = true;      end             % MMP, 09/25/2026
for j = find(~any(occ|zfill,1)),    zfill(1,j) = true;      end             % MMP, 09/25/2026
deg0 = struct('mult',0,'int',[0 0],'out',0,'in',0);                         % MMP, 09/25/2026

% % % First pass: each block's bases and the coefficient positions of its
% % % variables, and the total count, so ONE declaration covers them all.
spec = cell(M,N);       q = 0;
for i = 1:M
    for j = 1:N
%       if ~occ(i,j),   continue,   end                                     % MMP, 09/25/2026 (was)
        if ~occ(i,j) && ~zfill(i,j),    continue,   end                     % MMP, 09/25/2026
        dg = degs{i,j};     if zfill(i,j),  dg = deg0;  end                 % MMP, 09/25/2026
%       s = block_spec(sp_out{i},sp_in{j},d_out(i),d_in(j),degs{i,j},vars,dom); % MMP, 09/25/2026 (was)
        s = block_spec(sp_out{i},sp_in{j},d_out(i),d_in(j),dg,vars,dom);    % MMP, 09/25/2026
        if zfill(i,j)                                                       % MMP, 09/25/2026
            s.pos = cellfun(@(p) zeros(0,1),s.pos,'UniformOutput',false);   % MMP, 09/25/2026
        end                                                                 % MMP, 09/25/2026
        s.first = q+1;
        % s.pos is 3 x 3 x ... for two or more shared variables, so sum
        % over ALL of it; a plain sum() returns a row and q a vector.
        q = q + sum(reshape(cellfun(@numel,s.pos),[],1));
        spec{i,j} = s;
    end
end

Zd = cell(0,1);
if q>0
    [prog,dv] = lpidecvar(prog,[q,1]);
    Zd = reshape(dv.dvarname,[],1);
    if numel(Zd)~=q
        error("Internal error: declared %d decision variables, expected %d.",numel(Zd),q)
    end
end

% % % Second pass: the blocks, straight onto the one list. Variable k of the
% % % container sits at one coefficient position, so each gamma cell's B is
% % % one 'sparse' from triplets with a single nonzero per variable.
C = cell(M,N);
for i = 1:M
    for j = 1:N
%       if ~occ(i,j),   continue,   end                                     % MMP, 09/25/2026 (was)
        if isempty(spec{i,j}),  continue,   end     % a [] block            % MMP, 09/25/2026
        s = spec{i,j};
        A = cell(size(s.pos));      B = cell(size(s.pos));
        k0 = s.first;
        for g = 1:numel(s.pos)
            ng = numel(s.pos{g});
            A{g} = sparse(s.nC,1);
            B{g} = sparse(k0+(0:ng-1)',s.pos{g},1,q,s.nC);
            k0 = k0+ng;
        end
        params = struct('A',{A},'B',{B});
        C{i,j} = sdopvar(params,s.vars,Zd,s.ZL,s.ZR,s.dom,s.dims);
    end
end

% % % The container, with its metadata stated: a structurally zero block
% % % leaves nothing to derive a row's or column's space from.
so = false(M,nv);       si = false(N,nv);
for i = 1:M,    so(i,:) = ismember(vars,sp_out{i});    end
for j = 1:N,    si(j,:) = ismember(vars,sp_in{j});     end
meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
    'dim_out',d_out(:),'dim_in',d_in(:),'Zd',{Zd});
Pop = cdopvar(C,meta);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = block_spec(vo,vi,m,n,dg,vars,dom)
% Bases, canonical variable order and coefficient positions of one block,
% from output space vo (m components) to input space vi (n components).

vo = sort(reshape(vo,1,[]));    vi = sort(reshape(vi,1,[]));
S3 = intersect(vo,vi);          S2 = setdiff(vo,vi);    S1 = setdiff(vi,vo);
S3 = reshape(S3,1,[]);          S2 = reshape(S2,1,[]);  S1 = reshape(S1,1,[]);
n1 = numel(S1);     n2 = numel(S2);     n3 = numel(S3);
vout = [S2, S3];    vin = [S3, S1];                 % the class's canonical order

% A shared variable's left basis must hold both the multiplier degrees and
% the integral ones; the right basis only the integral ones, since a
% multiplier direction carries no dummy variable degree.
dl3 = max(dg.mult,dg.int(1));
ZL = [repmat({(0:dg.out)'},1,n2), repmat({(0:dl3)'},1,n3)];
ZR = [repmat({(0:dg.int(2))'},1,n3), repmat({(0:dg.in)'},1,n1)];
nL = cellfun(@numel,ZL);    nR = cellfun(@numel,ZR);
NL = prod([nL,1]);          NR = prod([nR,1]);
nrow = m*NL;                nC = nrow*n*NR;
sL = strides(nL);           sR = strides(nR);

ncell = 3^n3;
pos = cell([3*ones(1,n3),1,1]);
for g = 1:ncell
    gam = gamma_of(g,n3);
    % Admissible left and right monomials, per variable, as 0-based indices.
    iL = cell(1,n2+n3);
    for p = 1:n2,   iL{p} = 0:nL(p)-1;      end
    for t = 1:n3
        cap = dg.int(1);    if gam(t)==1,   cap = dg.mult;  end
        iL{n2+t} = find(ZL{n2+t}<=cap)'-1;
    end
    iR = cell(1,n3+n1);
    for t = 1:n3
        cap = dg.int(2);    if gam(t)==1,   cap = 0;        end
        iR{t} = find(ZR{t}<=cap)'-1;
    end
    for p = 1:n1,   iR{n3+p} = 0:nR(n3+p)-1;    end
    a = tensor_index(iL,sL);        b = tensor_index(iR,sR);
    % Coefficient rows are (component, left monomial) with the component
    % outer, columns (component, right monomial) likewise; vec is column
    % major - the layout 'canonicalize_multiplier' reads.
    r = reshape((0:m-1)'*NL + reshape(a,1,[]),[],1);
    c = reshape((0:n-1)'*NR + reshape(b,1,[]),[],1);
    pos{g} = reshape(r + nrow*reshape(c,1,[]),[],1) + 1;
end

[~,io] = ismember(vout,vars);   [~,ii] = ismember(vin,vars);
s = struct();
s.vars = struct('out',{vout},'in',{vin});
s.dom  = struct('out',dom(io,:),'in',dom(ii,:));
s.ZL = ZL;      s.ZR = ZR;      s.dims = [m,n];
s.nC = nC;      s.pos = pos;

end


function st = strides(nvec)
% Stride of each variable in kron(Z{1},...,Z{N}), first variable slowest.
st = ones(1,numel(nvec));
for k = numel(nvec)-1:-1:1,     st(k) = st(k+1)*nvec(k+1);   end
end


function lin = tensor_index(idx,st)
% Linear 0-based indices of the tensor product of per-variable index sets.
lin = 0;
for p = 1:numel(idx)
    lin = reshape(lin(:) + st(p)*reshape(idx{p},1,[]),[],1);
end
end


function gam = gamma_of(g,n3)
% The gamma multi-index of linear cell g, first direction fastest, as the
% class's ind2sub over [3 ... 3] gives it.
if n3==0,   gam = zeros(1,0);   return,     end
c = cell(1,n3);
[c{:}] = ind2sub([3*ones(1,n3),1],g);
gam = cell2mat(c(1:n3));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [so,si,dout,din] = parse_spaces(dims,spaces)
% Output and input space lists, and component counts.
if isstruct(spaces)
    if ~isscalar(spaces)
        error("'spaces' should be a single struct, not an array; struct() with "...
              +"a cell value returns an array, so assign the fields instead.")
    end
    if ~isfield(spaces,'out')
        error("A 'struct' spaces argument needs a field 'out'.")
    end
    so = norm_spaces(spaces.out);
    if isfield(spaces,'in'),    si = norm_spaces(spaces.in);
    else,                       si = so;
    end
else
    so = norm_spaces(spaces);   si = so;
end
M = numel(so);      N = numel(si);
if isstruct(dims)
    dout = reshape(dims.out,[],1);
    if isfield(dims,'in'),  din = reshape(dims.in,[],1);
    else,                   din = dout;
    end
else
    dout = reshape(dims,[],1);  din = dout;
end
if isscalar(dout),  dout = repmat(dout,M,1);    end
if isscalar(din),   din = repmat(din,N,1);      end
if numel(dout)~=M || numel(din)~=N
    error("Dimensions should have one entry per space: %d output and %d input.",M,N)
end
if any([dout;din]<1) || any([dout;din]~=round([dout;din]))
    error("Component counts should be positive integers.")
end
end


function sp = norm_spaces(sp)
% A plain cellstr is ONE space; otherwise a cell of cellstr, one per space.
if isempty(sp),     sp = {cell(1,0)};   return,     end
if iscellstr(sp) || ischar(sp) || isa(sp,'polynomial'),  sp = {sp};  end
sp = reshape(sp,1,[]);
for k = 1:numel(sp)
    sk = sp{k};
    if isa(sk,'polynomial'),            sk = sk.varname(:)';    end
    if ischar(sk),                      sk = {sk};              end
    if isnumeric(sk) && isempty(sk),    sk = cell(1,0);         end
    if ~iscellstr(sk)
        error("Space %d should be a 'cellstr' of variable names; an empty "...
              +"one is the finite-dimensional space R^q.",k)
    end
    sk = reshape(sk,1,[]);
    if numel(unique(sk))~=numel(sk)
        error("Space %d repeats a variable name.",k)
    end
    sp{k} = sk;
end
end


function dom = parse_dom(dom,vars)
% As 'copquadvar': nv x 2 in registry order, one 1 x 2 row for all, or a
% struct pairing names with intervals.
nv = numel(vars);
if isa(dom,'struct')
    if ~isfield(dom,'vars') || ~isfield(dom,'dom')
        error("A 'struct' domain should have fields 'vars' and 'dom'.")
    end
    if ~isscalar(dom)
        error("A 'struct' domain should be a single struct, not an array; "...
              +"struct() with a cell value returns an array.")
    end
    dvars = dom.vars;
    if isa(dvars,'polynomial'),  dvars = dvars.varname(:)';   end
    if ischar(dvars),            dvars = {dvars};             end
    dvars = reshape(dvars,1,[]);
    if size(dom.dom,1)~=numel(dvars) || size(dom.dom,2)~=2
        error("A 'struct' domain should pair each name in 'vars' with a row of 'dom'.")
    end
    [tf,loc] = ismember(vars,dvars);
    if ~all(tf)
        error("No domain was given for the variable '"+string(vars{find(~tf,1)})+"'.")
    end
    dom = dom.dom(loc,:);
end
if nv==0
    dom = zeros(0,2);
    return
end
if size(dom,2)~=2
    error("Domains should be specified as an nv x 2 array.")
end
if size(dom,1)==1 && nv~=1
    dom = repmat(dom,nv,1);
elseif size(dom,1)~=nv
    error("Domains should be an nv x 2 array for the %d registry variables.",nv)
end
if any(dom(:,2)<=dom(:,1))
    error("Each domain should satisfy dom(d,1) < dom(d,2).")
end
end


function degs = parse_deg(deg,M,N)
% One normalized specification per block: fields mult, int (1 x 2), out, in.
if iscell(deg)
    if ~isequal(size(deg),[M,N])
        error("A cell of degree specifications should be %d x %d.",M,N)
    end
    degs = cellfun(@norm_deg,deg,'UniformOutput',false);
else
    d = norm_deg(deg);
    degs = repmat({d},M,N);
end
end


function d = norm_deg(x)
% 'get_lpivar_degs' returns its degrees SPARSE, and the executives pass them
% straight on; ranges built from a sparse scalar are sparse, so make it full.
if isnumeric(x),    x = full(x);    end
if isnumeric(x) && isscalar(x)
    d = struct('mult',x,'int',[x,x],'out',x,'in',x);
elseif isnumeric(x) && numel(x)==3
    d = struct('mult',x(1),'int',[x(3),x(2)],'out',x(1),'in',x(1));
elseif isstruct(x) && isscalar(x)
    d = struct('mult',1,'int',[1,1],'out',1,'in',1);
    f = fieldnames(x);
    for k = 1:numel(f)
        if ~isfield(d,f{k})
            error("Unknown degree field '%s'; use mult, int, out or in.",f{k})
        end
        d.(f{k}) = x.(f{k});
    end
    if isscalar(d.int),     d.int = [d.int,d.int];      end
else
    error("A degree should be a scalar, a 1x3 array [d1 d2 d3], or a struct "...
          +"with fields mult, int, out, in.")
end
v = [d.mult, d.int(:)', d.out, d.in];
if numel(d.int)~=2 || any(v<0) || any(v~=round(v))
    error("Degrees should be nonnegative integers, and 'int' a 1 x 2 [left right].")
end
end
