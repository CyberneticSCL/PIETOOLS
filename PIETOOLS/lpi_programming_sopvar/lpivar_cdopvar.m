function [prog,Pop] = lpivar_cdopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP] = LPIVAR_CDOPVAR(PROG,DIMS,SPACES,DOM,DEG,OPTIONS) declares in
% the LPI program PROG the decision operator
%
%   Pop: L_2^{p_1}[s^1] x ... x L_2^{p_N}[s^N]
%                   -> L_2^{q_1}[t^1] x ... x L_2^{q_M}[t^M],
%   (Pop x)_i = sum_j Pop.C{i,j} x_j,
%
% with free kernels. Block (i,j) maps input space j to output space i;
% with S3 the variables the two spaces share, S2 those of the output space
% only and S1 those of the input space only, its kernel in the parameter
% cell gam in {1,2,3}^|S3| is
%
%   K_gam(t,s') = (I_{q_i} (x) ZL(t))' C_gam (I_{p_j} (x) ZR(s')),
%
% with t the variables [S2,S3], s' the variables [S3,S1], ZL and ZR the
% Kronecker products of the monomial vectors of each variable, and, for
% the k-th variable of S3, gam_k = 1 a multiplier delta(t_k - s'_k), 2 and
% 3 the integrals over s'_k < t_k and s'_k > t_k. Every entry of C_gam
% whose monomials are within the caps of DEG, and of degree 0 in s'_k
% where gam_k = 1, is a separate decision variable; every other entry of
% C_gam is zero.
%
% [PROG,POP] = LPIVAR_CDOPVAR(PROG,X,DEG,OPTIONS) declares Pop with the
% output spaces, the input spaces and the domain of the 'copvar' or
% 'cdopvar' object X, so that Pop has the shape of X.
%
% INPUT
% - prog:    'struct', the LPI program, with the field decvartable;
% - dims:    'struct' with fields out (M x 1 'double', the counts q_i) and,
%            optionally, in (N x 1 'double', the counts p_j; default out);
%            or a 'double' array, the counts q_i, also taken as the p_j. A
%            scalar count applies to every space. Counts are positive
%            integers. In the 'struct' form an empty out or in is no space
%            on that side (an empty channel of a PIE), and SPACES of that
%            side is not read; an empty 'double' array is no space on
%            either side;
% - spaces:  'struct' with fields out (1 x M cell, entry i the 'cellstr' of
%            the variable names of t^i) and, optionally, in (1 x N cell,
%            the same for s^j; default out); or a 1 x M cell, the output
%            spaces, also taken as the input spaces. An empty entry is the
%            finite-dimensional space R^{q_i} (R^{p_j}); a 'cellstr', char
%            or 'polynomial' in place of the cell is one space. No space
%            repeats a name, and no name ends in '_int' or '_dum';
% - dom:     nv x 2 'double', row k the interval [a,b], a < b, of the k-th
%            name in the sorted list of the nv names of SPACES; one 1 x 2
%            row for every variable; or a 'struct' with fields vars
%            ('cellstr', every name of SPACES among them) and dom
%            ('double', one row [a,b] per name of vars). Optional: omitted,
%            as LPIVAR_CDOPVAR(PROG,DIMS,SPACES,DEG,OPTIONS), the domain
%            of PROG, by name;
% - X:       'copvar' or 'cdopvar' object, in place of DIMS, SPACES and
%            DOM (second form above);
% - deg:     (optional) the degree caps, nonnegative integers (default: per
%            block, from lpi.var.deg in the program settings, the left
%            caps mult, int(1) and out at the entry of the dimension of the
%            output space, the right caps int(2) and in at that of the
%            input space; 1 for a program without settings):
%            a 'struct' with any of the fields
%              mult  'double', the degree in t_k where gam_k = 1;
%              int   1 x 2 'double', the degrees in t_k and in s'_k where
%                    gam_k = 2 or 3 (a scalar for both);
%              out   'double', the degree in each variable of S2;
%              in    'double', the degree in each variable of S1;
%            a field not given being 1; or a 'double' d, for mult = out =
%            in = d and int = [d d]; or a 1 x 3 'double' [d1 d2 d3], for
%            mult = out = in = d1 and int = [d3 d2]; or an M x N cell of
%            these, entry (i,j) for block (i,j);
% - options: (optional) 'struct' with the optional field occ (M x N
%            'logical', default true(M,N)): occ(i,j) false makes block
%            (i,j) structurally zero ([], no decision variables); a row or
%            column with no true entry gets, in its first column or row, a
%            block of degree 0 with zero coefficients and no decision
%            variables.
%
% OUTPUT
% - prog:    'struct', PROG with the decision variables of POP appended to
%            prog.decvartable;
% - Pop:     'cdopvar' object of M x N blocks, block (i,j) an 'sdopvar'
%            object with output variables [S2,S3] and input variables
%            [S3,S1] (each set sorted) and the exponents 0..out (S2) and
%            0..max(mult,int(1)) (S3) in ZL, 0..int(2) (S3) and 0..in (S1)
%            in ZR; [] where occ is false, except the blocks of degree 0
%            above. Pop.Zd, the 'cellstr' column of the decision variable
%            names, is shared by every block. With no space on a side, POP
%            is M x 0 or 0 x N (0 x 0 if both), with no block and no
%            decision variable, and PROG is returned unchanged.
%
% See also POSCOPVAR, LPI_EQ_CDOPVAR, CDOPVAR, LPIPROGRAM.
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

% Defaults of the optional inputs; PROG must be an LPI program
% DOM omitted: the 4th input, if any, is DEG (a scalar, 1 x 3, a 'struct'
% without field dom, or a 'cell', never n x 2), and the domain is that of
% PROG by name, per element, since the varname of a polynomial vector is
% sorted while the rows of prog.dom follow the program order
nin = nargin;
% Operator form LPIVAR_CDOPVAR(PROG,X,DEG,OPTIONS): the output and input
% spaces of X and its domain stand for DIMS, SPACES, DOM; DEG and OPTIONS
% shift by two. Any other DIMS: the original form
if isa(dims,'copvar') || isa(dims,'cdopvar')
    narginchk(2,4);
    if nin>=4,  options = dom;  end
    if nin>=3,  deg = spaces;   end
    X = dims;   dims = struct();    spaces = struct();
    [spaces.out,dims.out] = copvar_space_list(X,'out');
    [spaces.in,dims.in] = copvar_space_list(X,'in');
    dom = struct('vars',{X.vars},'dom',X.dom);  nin = nin+2;
end
if nin<4 || ~((isnumeric(dom) && (isempty(dom) || size(dom,2)==2)) ...
              || (isstruct(dom) && isfield(dom,'dom')))
    if nin>=5,  options = deg;  end
    if nin>=4,  deg = dom;      end
    nin = nin+1;
    v = cell(1,size(prog.dom,1));
    for k = 1:numel(v),     v{k} = prog.vartable(k).varname{1};     end
    dom = struct('vars',{v},'dom',prog.dom);
end
defdeg = nin<5 || isempty(deg);     % from the settings, after the parse
if nin<6 || isempty(options),       options = struct(); end
if ~isa(prog,'struct') || ~isfield(prog,'decvartable')
    error("'prog' should be an LPI program structure; use 'lpiprogram'.")
end

% An empty count list on a side: no space there, so no block and no
% decision variable; PROG is unchanged
if (~isstruct(dims) && isempty(dims)) || ...
        (isstruct(dims) && ((isfield(dims,'out') && isempty(dims.out)) || ...
                            (isfield(dims,'in') && isempty(dims.in))))
    Pop = empty_side(dims,spaces,dom);
    return
end

% % % Spaces and dimensions

% % % Registry, sorted, as the containers keep it: an sdopvar indexes its
% gamma cells over S3 in sorted order, and with a sorted registry every
% space and shared set is already sorted. Parse, in this order: spaces,
% dims, reserved '_int'/'_dum' names, domain.
[meta,sp_out,sp_in] = parse_copvar_spaces(dims,spaces,dom);
M = numel(sp_out);      N = numel(sp_in);
% DEG omitted: lpi.var.deg of the program settings, per block (i,j): the
% left caps (mult, int(1), out) from the entry of the dimension of output
% space i, the right caps (int(2), in) from that of input space j, so a
% block depends on its own two spaces only. R^n reads entry 1, which then
% caps no variable. 1 for a program without settings.
if defdeg
    deg = 1;
    if isfield(prog,'settings') && isfield(prog.settings,'lpi') && isfield(prog.settings.lpi,'var')
        zd = prog.settings.lpi.var.deg;
        e = @(a) norm_deg(zd{min(max(a,1),end)});
        dO = arrayfun(e,sum(meta.space_out,2),'UniformOutput',false);
        dI = arrayfun(e,sum(meta.space_in,2),'UniformOutput',false);
        deg = cell(M,N);
        for i = 1:M
            for j = 1:N
                deg{i,j} = struct('mult',dO{i}.mult,'int',[dO{i}.int(1),dI{j}.int(2)], ...
                                  'out',dO{i}.out,'in',dI{j}.in);
            end
        end
    end
end
vars = meta.vars;       dom = meta.dom;
d_out = meta.dim_out;   d_in = meta.dim_in;

% % % Occupancy and degrees
occ = true(M,N);
if isfield(options,'occ') && ~isempty(options.occ)
    occ = logical(options.occ);
    if ~isequal(size(occ),[M,N])
        error("options.occ should be %d x %d, one entry per block.",M,N)
    end
end
% The degree caps of every block in normalized form, fields mult, int
% (1 x 2), out, in: one specification per block from an M x N cell, or one
% specification for every block
if iscell(deg)
    if ~isequal(size(deg),[M,N])
        error("A cell of degree specifications should be %d x %d.",M,N)
    end
    degs = cellfun(@norm_deg,deg,'UniformOutput',false);
else
    d = norm_deg(deg);
    degs = repmat({d},M,N);
end
% A row or column with no block cannot have its space read off one, so it
% gets one explicit zero block: degree 0, no coefficient positions, hence
% no variables. Rows first, then columns.
zfill = false(M,N);
for i = find(~any(occ,2))',         zfill(i,1) = true;      end
for j = find(~any(occ|zfill,1)),    zfill(1,j) = true;      end
deg0 = struct('mult',0,'int',[0 0],'out',0,'in',0);

% % % First pass: the bases of each block and the coefficient positions of
% % % its variables, and the total count, so ONE declaration covers them all.
spec = cell(M,N);       q = 0;
for i = 1:M
    for j = 1:N
        if ~occ(i,j) && ~zfill(i,j),    continue,   end
        dg = degs{i,j};     if zfill(i,j),  dg = deg0;  end
        s = block_spec(sp_out{i},sp_in{j},d_out(i),d_in(j),dg,vars,dom);
        if zfill(i,j)
            s.pos = cellfun(@(p) zeros(0,1),s.pos,'UniformOutput',false);
        end
        % Index of the first decision variable of the block on the one list
        s.first = q+1;
        % s.pos is 3 x 3 x ... for two or more shared variables, so sum
        % over ALL of it; a plain sum() returns a row and q a vector.
        q = q + sum(reshape(cellfun(@numel,s.pos),[],1));
        spec{i,j} = s;
    end
end

% One declaration of the q decision variables; Zd their names
Zd = cell(0,1);
if q>0
    [prog,dv] = lpidecvar(prog,[q,1]);
    Zd = reshape(dv.dvarname,[],1);
    if numel(Zd)~=q
        error("Internal error: declared %d decision variables, expected %d.",numel(Zd),q)
    end
end

% % % Second pass: the blocks, straight onto the one list. Variable k of the
% % % container sits at one coefficient position, so the B of each gamma
% % % cell is one 'sparse' from triplets with a single nonzero per variable.
C = cell(M,N);
for i = 1:M
    for j = 1:N
        if isempty(spec{i,j}),  continue,   end     % a [] block
        s = spec{i,j};
        % Gamma cell g: a zero constant part, and in B a 1 in row k0+r-1
        % and column s.pos{g}(r): decision variable k0+r-1 is the
        % coefficient at that position of vec(C_g)
        A = cell(size(s.pos));      B = cell(size(s.pos));
        k0 = s.first;
        for g = 1:numel(s.pos)
            ng = numel(s.pos{g});
            A{g} = sparse(s.nC,1);
            B{g} = sparse(k0+(0:ng-1)',s.pos{g},1,q,s.nC);
            k0 = k0+ng;
        end
        % The block, on the decision variable list Zd of the container
        params = struct('A',{A},'B',{B});
        C{i,j} = sdopvar(params,s.vars,Zd,s.ZL,s.ZR,s.dom,s.dims);
    end
end

% % % The container, with its metadata stated: a structurally zero block
% % % leaves nothing to derive the space of a row or column from.
% meta holds these fields already, in this order; Zd is added last. The
% assignment shares Zd, it does not copy q names.
meta.Zd = Zd;
Pop = cdopvar(C,meta);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = block_spec(vo,vi,m,n,dg,vars,dom)
% S = BLOCK_SPEC(VO,VI,M,N,DG,VARS,DOM) returns the variables, domains,
% monomial bases and decision coefficient positions of the block from the
% input space VI (N components) to the output space VO (M components).
% INPUT  vo, vi: cellstr, the variable names of the output and input space;
%        m, n:   double, their component counts;
%        dg:     'struct' of degree caps, fields mult, int (1 x 2), out, in;
%        vars:   1 x nv cellstr, the sorted registry;
%        dom:    nv x 2 double, the domains in registry order.
% OUTPUT s:      'struct' with fields vars and dom ('struct's with fields
%                out and in), ZL and ZR (cells of exponent columns), dims
%                ([m,n]), nC (double, the entries of a coefficient matrix)
%                and pos (3^n3 cells, 3 x ... x 3; cell g the 1-based
%                positions in vec(C_g) of the admissible coefficients).

% Shared (S3), output-only (S2) and input-only (S1) variables, each sorted
vo = sort(reshape(vo,1,[]));    vi = sort(reshape(vi,1,[]));
S3 = intersect(vo,vi);          S2 = setdiff(vo,vi);    S1 = setdiff(vi,vo);
S3 = reshape(S3,1,[]);          S2 = reshape(S2,1,[]);  S1 = reshape(S1,1,[]);
n1 = numel(S1);     n2 = numel(S2);     n3 = numel(S3);
vout = [S2, S3];    vin = [S3, S1];                 % the canonical order of the class

% The left basis of a shared variable must hold both the multiplier
% degrees and the integral ones; the right basis only the integral ones,
% since a multiplier direction carries no dummy variable degree.
dl3 = max(dg.mult,dg.int(1));
ZL = [repmat({(0:dg.out)'},1,n2), repmat({(0:dl3)'},1,n3)];
ZR = [repmat({(0:dg.int(2))'},1,n3), repmat({(0:dg.in)'},1,n1)];
% Basis sizes, the coefficient matrix (nrow x n*NR, nC entries), and the
% Kronecker strides of the two bases
nL = cellfun(@numel,ZL);    nR = cellfun(@numel,ZR);
NL = prod([nL,1]);          NR = prod([nR,1]);
nrow = m*NL;                nC = nrow*n*NR;
sL = kron_strides(nL);      sR = kron_strides(nR);

% Positions of the admissible coefficients, one gamma cell at a time
ncell = 3^n3;
pos = cell([3*ones(1,n3),1,1]);
for g = 1:ncell
    gam = gamma_of_cell(g,n3);      % direction 1 fastest, as the class
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
    % 0-based Kronecker positions of the product set, variable 1 fastest in
    % the enumeration: exact integer arithmetic.
    a = multiindex_grid(iL)*sL(:);  b = multiindex_grid(iR)*sR(:);
    % Coefficient rows are (component, left monomial) with the component
    % outer, columns (component, right monomial) likewise; vec is column
    % major.
    r = reshape((0:m-1)'*NL + reshape(a,1,[]),[],1);
    c = reshape((0:n-1)'*NR + reshape(b,1,[]),[],1);
    pos{g} = reshape(r + nrow*reshape(c,1,[]),[],1) + 1;
end

% The variables and domains of the block on each side, bases, size,
% positions
[~,io] = ismember(vout,vars);   [~,ii] = ismember(vin,vars);
s = struct();
s.vars = struct('out',{vout},'in',{vin});
s.dom  = struct('out',dom(io,:),'in',dom(ii,:));
s.ZL = ZL;      s.ZR = ZR;      s.dims = [m,n];
s.nC = nC;      s.pos = pos;

end


function d = norm_deg(x)
% D = NORM_DEG(X) returns one degree specification as a 'struct' with the
% fields mult, int (1 x 2), out and in, each a nonnegative integer.
% INPUT  x: 'double' d (every cap d), 1 x 3 'double' [d1 d2 d3] (mult =
%           out = in = d1, int = [d3 d2]), or 'struct' with any of the
%           fields mult, int, out, in (the others 1).
% OUTPUT d: 'struct'.
% Degrees may come SPARSE; ranges built from a sparse scalar are sparse, so
% make it full.
if isnumeric(x),    x = full(x);    end
% A scalar, a 1 x 3 array, or a struct over the defaults
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
% Every cap a nonnegative integer, and 'int' of two entries
v = [d.mult, d.int(:)', d.out, d.in];
if numel(d.int)~=2 || any(v<0) || any(v~=round(v))
    error("Degrees should be nonnegative integers, and 'int' a 1 x 2 [left right].")
end
end


function Pop = empty_side(dims,spaces,dom)
% POP = EMPTY_SIDE(DIMS,SPACES,DOM) returns the container with no space on
% the side whose count list in DIMS is empty and the spaces of the other
% side, with no block and no decision variable.
% INPUT  dims:   'struct' with fields out and (optionally) in, at least
%                one of them empty (in absent: the same as out), or an
%                empty 'double' (no space on either side);
%        spaces: 'struct' with fields out and in, as LPIVAR_CDOPVAR; the
%                field of an empty side is not read;
%        dom:    as LPIVAR_CDOPVAR.
% OUTPUT Pop:    'cdopvar' object, M x 0, 0 x N or 0 x 0, Zd empty.
% Which sides are empty; an empty 'double' DIMS is both
if ~isstruct(dims),     dims = struct('out',[],'in',[]);    end
eo = isempty(dims.out);
ei = eo;    if isfield(dims,'in'),  ei = isempty(dims.in);  end
if eo && ei
    % Neither side has a space
    meta = struct('vars',{cell(1,0)},'dom',zeros(0,2),'space_out',false(0,0),...
                  'space_in',false(0,0),'dim_out',zeros(0,1),'dim_in',zeros(0,1));
else
    % The spaces of the other side, parsed as one list (output = input),
    % then the empty side stated with no space
    if eo,  meta = parse_copvar_spaces(dims.in,spaces.in,dom);
    else,   meta = parse_copvar_spaces(dims.out,spaces.out,dom);
    end
    nv = numel(meta.vars);
    if eo
        meta.space_in = meta.space_out;     meta.dim_in = meta.dim_out;
        meta.space_out = false(0,nv);       meta.dim_out = zeros(0,1);
    else
        meta.space_in = false(0,nv);        meta.dim_in = zeros(0,1);
    end
end
% No block, no decision variable
meta.Zd = cell(0,1);
Pop = cdopvar(cell(numel(meta.dim_out),numel(meta.dim_in)),meta);
end
