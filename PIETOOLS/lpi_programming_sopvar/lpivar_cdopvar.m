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
% INPUT
% - prog:    'struct', the LPI program (LPIPROGRAM_SOP), with the field
%            decvartable;
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
%            ('double', one row [a,b] per name of vars);
% - deg:     (optional) the degree caps, nonnegative integers (default 1):
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
% MMP, 09/30/2026: Shared helpers of sopvar/misc/conventions replace local
%                  copies of conventions that three to four files each
%                  wrote: 'parse_copvar_spaces' reads spaces, dimensions,
%                  sorted registry, reserved names and domains (its bodies
%                  are this file's parse_spaces, norm_spaces and parse_dom,
%                  which also lived verbatim in 'spaces2meta_sop');
%                  'kron_strides', 'gamma_of_cell' and 'multiindex_grid'
%                  give the Kronecker strides, the gamma of a parameter cell
%                  and the tensor indices (the cells: Sec. 1 and 4 of
%                  sopvar_implementation_notes.pdf).
%                  Locals strides, tensor_index, gamma_of, parse_spaces,
%                  norm_spaces and parse_dom deleted. Same program bit for
%                  bit (tr_ab, 91 SDP leaves over w1-w3). One message
%                  changes: a non-cell space list, e.g. 5, now raises the
%                  parser's "Spaces should be specified as a cell of
%                  'cellstr' objects." instead of MATLAB's brace-indexing
%                  error. NOTES line on 'copquadvar' replaced.
% MMP, 09/30/2026: Help for 'deg' states that 'int' and 'mult' are kernel
%                  degrees here but basis-operator caps in 'copquadvar' and
%                  'sopquadvar', and that a struct with only those two
%                  fields, 'mult' scalar, passes both (measured at nv = 1
%                  and 2), so the mix-up raises no error. Help only; the
%                  field names are unchanged.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body and headers
%                  of the local functions; code unchanged. The deleted
%                  header text was the comparisons with lpivar, lpivar_2d,
%                  copquadvar and poscopvar (among them the 'deg' remark of
%                  the third 09/30/2026 entry), the variable count and span
%                  against lpivar with a test reference, the note on the
%                  shared parser of the second 09/30/2026 entry, and a cost
%                  note.
% MMP, 10/28/2026: An empty count list on one side of a 'struct' DIMS is no
%                  space on that side: an M x 0 or 0 x N container with no
%                  decision variable (new local function empty_side).
%                  Before, the parser read the empty space list as one R^q
%                  space and the counts as missing, so the container
%                  executives failed on a PIE with an empty channel (the
%                  observer gain from no measurement, the controller of no
%                  input); maintainer: empty channels should not raise an
%                  error. Help text: DIMS and POP.

% Defaults of the optional inputs; PROG must be an LPI program
if nargin<5 || isempty(deg),        deg = 1;            end
if nargin<6 || isempty(options),    options = struct(); end
if ~isa(prog,'struct') || ~isfield(prog,'decvartable')
    error("'prog' should be an LPI program structure; use 'lpiprogram'.")
end

% An empty count list on a side: no space there, so no block and no         % MMP, 10/28/2026
% decision variable; PROG is unchanged                                      % MMP, 10/28/2026
if (~isstruct(dims) && isempty(dims)) || ...
        (isstruct(dims) && ((isfield(dims,'out') && isempty(dims.out)) || ...
                            (isfield(dims,'in') && isempty(dims.in))))      % MMP, 10/28/2026
    Pop = empty_side(dims,spaces,dom);                                      % MMP, 10/28/2026
    return                                                                  % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026

% % % Spaces and dimensions
% [sp_out,sp_in,d_out,d_in] = parse_spaces(dims,spaces);                    % MMP, 09/30/2026 (was)
% M = numel(sp_out);      N = numel(sp_in);                                 % MMP, 09/30/2026 (was)

% % % Registry, sorted, as the containers and 'poscopvar' keep it: an
% sdopvar indexes its gamma cells over S3 in sorted order, and with a sorted
% registry every space and shared set is already sorted.
% vars = reshape(unique([sp_out{:}, sp_in{:}]),1,[]);     nv = numel(vars); % MMP, 09/30/2026 (was)
% is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));       % MMP, 09/30/2026 (was)
% if any(is_reserved)                                                       % MMP, 09/30/2026 (was)
%     error("Spatial variable names may not end in '_int' or '_dum'; "...
%           +"'"+string(vars{find(is_reserved,1)})+"' does.")               % MMP, 09/30/2026 (was)
% end                                                                       % MMP, 09/30/2026 (was)
% dom = parse_dom(dom,vars);                                                % MMP, 09/30/2026 (was)
% One shared parser for these inputs; same values and errors, in the same   % MMP, 09/30/2026
% order: spaces, dims, reserved '_int'/'_dum' names, domain.                % MMP, 09/30/2026
[meta,sp_out,sp_in] = parse_copvar_spaces(dims,spaces,dom);                 % MMP, 09/30/2026
M = numel(sp_out);      N = numel(sp_in);                                   % MMP, 09/30/2026
vars = meta.vars;       dom = meta.dom;                                     % MMP, 09/30/2026
d_out = meta.dim_out;   d_in = meta.dim_in;                                 % MMP, 09/30/2026

% % % Occupancy and degrees
occ = true(M,N);
if isfield(options,'occ') && ~isempty(options.occ)
    occ = logical(options.occ);
    if ~isequal(size(occ),[M,N])
        error("options.occ should be %d x %d, one entry per block.",M,N)
    end
end
% Degree caps of every block
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
        % Index of the block's first decision variable on the one list
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
% % % container sits at one coefficient position, so each gamma cell's B is
% % % one 'sparse' from triplets with a single nonzero per variable.
C = cell(M,N);
for i = 1:M
    for j = 1:N
%       if ~occ(i,j),   continue,   end                                     % MMP, 09/25/2026 (was)
        if isempty(spec{i,j}),  continue,   end     % a [] block            % MMP, 09/25/2026
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
        % The block, on the container's decision variable list Zd
        params = struct('A',{A},'B',{B});
        C{i,j} = sdopvar(params,s.vars,Zd,s.ZL,s.ZR,s.dom,s.dims);
    end
end

% % % The container, with its metadata stated: a structurally zero block
% % % leaves nothing to derive a row's or column's space from.
% so = false(M,nv);       si = false(N,nv);                                 % MMP, 09/30/2026 (was)
% for i = 1:M,    so(i,:) = ismember(vars,sp_out{i});    end                % MMP, 09/30/2026 (was)
% for j = 1:N,    si(j,:) = ismember(vars,sp_in{j});     end                % MMP, 09/30/2026 (was)
% meta = struct('vars',{vars},'dom',dom,'space_out',so,'space_in',si,...
%     'dim_out',d_out(:),'dim_in',d_in(:),'Zd',{Zd});                       % MMP, 09/30/2026 (was)
% The parser's meta holds these fields already, in this order; Zd is added  % MMP, 09/30/2026
% last, as before. The assignment shares Zd, it does not copy q names.      % MMP, 09/30/2026
meta.Zd = Zd;                                                               % MMP, 09/30/2026
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
vout = [S2, S3];    vin = [S3, S1];                 % the class's canonical order

% A shared variable's left basis must hold both the multiplier degrees and
% the integral ones; the right basis only the integral ones, since a
% multiplier direction carries no dummy variable degree.
dl3 = max(dg.mult,dg.int(1));
ZL = [repmat({(0:dg.out)'},1,n2), repmat({(0:dl3)'},1,n3)];
ZR = [repmat({(0:dg.int(2))'},1,n3), repmat({(0:dg.in)'},1,n1)];
% Basis sizes, the coefficient matrix (nrow x n*NR, nC entries), and the
% Kronecker strides of the two bases
nL = cellfun(@numel,ZL);    nR = cellfun(@numel,ZR);
NL = prod([nL,1]);          NR = prod([nR,1]);
nrow = m*NL;                nC = nrow*n*NR;
% sL = strides(nL);           sR = strides(nR);                             % MMP, 09/30/2026 (was)
sL = kron_strides(nL);      sR = kron_strides(nR);                          % MMP, 09/30/2026

% Positions of the admissible coefficients, one gamma cell at a time
ncell = 3^n3;
pos = cell([3*ones(1,n3),1,1]);
for g = 1:ncell
%   gam = gamma_of(g,n3);                                                   % MMP, 09/30/2026 (was)
    gam = gamma_of_cell(g,n3);      % direction 1 fastest, as the class     % MMP, 09/30/2026
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
%   a = tensor_index(iL,sL);        b = tensor_index(iR,sR);                % MMP, 09/30/2026 (was)
    % 0-based Kronecker positions of the product set, variable 1 fastest in % MMP, 09/30/2026
    % the enumeration (as tensor_index): exact integer arithmetic.          % MMP, 09/30/2026
    a = multiindex_grid(iL)*sL(:);  b = multiindex_grid(iR)*sR(:);          % MMP, 09/30/2026
    % Coefficient rows are (component, left monomial) with the component
    % outer, columns (component, right monomial) likewise; vec is column
    % major - the layout 'canonicalize_multiplier' reads.
    r = reshape((0:m-1)'*NL + reshape(a,1,[]),[],1);
    c = reshape((0:n-1)'*NR + reshape(b,1,[]),[],1);
    pos{g} = reshape(r + nrow*reshape(c,1,[]),[],1) + 1;
end

% The block's variables and domains on each side, bases, size, positions
[~,io] = ismember(vout,vars);   [~,ii] = ismember(vin,vars);
s = struct();
s.vars = struct('out',{vout},'in',{vin});
s.dom  = struct('out',dom(io,:),'in',dom(ii,:));
s.ZL = ZL;      s.ZR = ZR;      s.dims = [m,n];
s.nC = nC;      s.pos = pos;

end


% MMP, 09/30/2026: local functions strides, tensor_index and gamma_of (Initial
% coding 09/25/2026, no line stamps) deleted, replaced by the shared
% kron_strides, multiindex_grid and gamma_of_cell (sopvar/misc/conventions),
% which reproduce them exactly (test_index_maps (B) holds verbatim copies).
% Local functions parse_spaces, norm_spaces and parse_dom (Initial coding
% 09/25/2026) deleted, moved verbatim into parse_copvar_spaces.


function degs = parse_deg(deg,M,N)
% DEGS = PARSE_DEG(DEG,M,N) returns the degree caps of each of the M x N
% blocks in normalized form.
% INPUT  deg:  'double' scalar or 1 x 3 array, 'struct', or M x N cell of
%              these (one per block);
%        M, N: double, the numbers of output and input spaces.
% OUTPUT degs: M x N cell of 'struct', fields mult, int (1 x 2), out, in.
if iscell(deg)
    % One specification per block
    if ~isequal(size(deg),[M,N])
        error("A cell of degree specifications should be %d x %d.",M,N)
    end
    degs = cellfun(@norm_deg,deg,'UniformOutput',false);
else
    % One specification for every block
    d = norm_deg(deg);
    degs = repmat({d},M,N);
end
end


function d = norm_deg(x)
% D = NORM_DEG(X) returns one degree specification as a 'struct' with the
% fields mult, int (1 x 2), out and in, each a nonnegative integer.
% INPUT  x: 'double' d (every cap d), 1 x 3 'double' [d1 d2 d3] (mult =
%           out = in = d1, int = [d3 d2]), or 'struct' with any of the
%           fields mult, int, out, in (the others 1).
% OUTPUT d: 'struct'.
% 'get_lpivar_degs' returns its degrees SPARSE, and the executives pass them
% straight on; ranges built from a sparse scalar are sparse, so make it full.
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


function Pop = empty_side(dims,spaces,dom)                                  % MMP, 10/28/2026
% POP = EMPTY_SIDE(DIMS,SPACES,DOM) returns the container with no space on  % MMP, 10/28/2026
% the side whose count list in DIMS is empty and the spaces of the other    % MMP, 10/28/2026
% side, with no block and no decision variable.                             % MMP, 10/28/2026
% INPUT  dims:   'struct' with fields out and (optionally) in, at least     % MMP, 10/28/2026
%                one of them empty (in absent: the same as out), or an      % MMP, 10/28/2026
%                empty 'double' (no space on either side);                  % MMP, 10/28/2026
%        spaces: 'struct' with fields out and in, as LPIVAR_CDOPVAR; the    % MMP, 10/28/2026
%                field of an empty side is not read;                        % MMP, 10/28/2026
%        dom:    as LPIVAR_CDOPVAR.                                         % MMP, 10/28/2026
% OUTPUT Pop:    'cdopvar' object, M x 0, 0 x N or 0 x 0, Zd empty.         % MMP, 10/28/2026
% Which sides are empty; an empty 'double' DIMS is both                     % MMP, 10/28/2026
if ~isstruct(dims),     dims = struct('out',[],'in',[]);    end             % MMP, 10/28/2026
eo = isempty(dims.out);                                                     % MMP, 10/28/2026
ei = eo;    if isfield(dims,'in'),  ei = isempty(dims.in);  end             % MMP, 10/28/2026
if eo && ei                                                                 % MMP, 10/28/2026
    % Neither side has a space                                              % MMP, 10/28/2026
    meta = struct('vars',{cell(1,0)},'dom',zeros(0,2),'space_out',false(0,0),...
                  'space_in',false(0,0),'dim_out',zeros(0,1),'dim_in',zeros(0,1)); % MMP, 10/28/2026
else                                                                        % MMP, 10/28/2026
    % The spaces of the other side, parsed as one list (output = input),    % MMP, 10/28/2026
    % then the empty side stated with no space                              % MMP, 10/28/2026
    if eo,  meta = parse_copvar_spaces(dims.in,spaces.in,dom);              % MMP, 10/28/2026
    else,   meta = parse_copvar_spaces(dims.out,spaces.out,dom);            % MMP, 10/28/2026
    end                                                                     % MMP, 10/28/2026
    nv = numel(meta.vars);                                                  % MMP, 10/28/2026
    if eo                                                                   % MMP, 10/28/2026
        meta.space_in = meta.space_out;     meta.dim_in = meta.dim_out;     % MMP, 10/28/2026
        meta.space_out = false(0,nv);       meta.dim_out = zeros(0,1);      % MMP, 10/28/2026
    else                                                                    % MMP, 10/28/2026
        meta.space_in = false(0,nv);        meta.dim_in = zeros(0,1);       % MMP, 10/28/2026
    end                                                                     % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
% No block, no decision variable                                            % MMP, 10/28/2026
meta.Zd = cell(0,1);                                                        % MMP, 10/28/2026
Pop = cdopvar(cell(numel(meta.dim_out),numel(meta.dim_in)),meta);           % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026
