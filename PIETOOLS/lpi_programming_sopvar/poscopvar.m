function [prog,Pop,Qcell,basis_list] = poscopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,BASIS_LIST] = POSCOPVAR(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares in the LPI program PROG a positive semidefinite, self-adjoint
% 'cdopvar' decision operator on a concatenation of spaces, parameterized
% as
%
%       Pop: X -> X,     X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M],
%
%       Pop = sum_{b,b'} Z_b'*(g*Q_{b,b'})*Z_{b'},     Q >= 0,
%
% with ONE Gram matrix Q of decision variables over every pair (b,b') of
% basis operators of all spaces, so that Pop is positive as a whole, and
% g(theta) >= 0 the weight of OPTIONS.psatz. An empty s^k is the space
% R^{m_k}. Basis operator b of space k carries a multi-index alpha over the
% variables of s^k and maps L_2^{m_k}[s^k] to L_2[theta], theta over every
% variable of the registry:
%
%       (Z_b x)(theta) = int I_alpha(theta-s) ...
%                            (I_{m_k} kron Z^alpha(theta,s)) x(s) ds,
%
% I_alpha(r) = prod_d I_{alpha_d}(r_d) over the variables of s^k, alpha_d =
% 1 the multiplier (delta), 2 the lower integral (s_d <= theta_d), 3 the
% upper integral (s_d >= theta_d), 4 the full integral; Z^alpha(theta,s)
% is the vector of monomials of the degrees DEG, and on R^{m_k},
% Z_b x = (I_{m_k} kron Z^alpha(theta)) x. Pop is also appended to the
% cell prog.sop_pos, the record of the positive operators of PROG. The
% position map of each distinct shape is cached for the MATLAB session
% (about 16 B per coefficient nonzero).
%
% [PROG,POP,QCELL,BASIS_LIST] = POSCOPVAR(PROG,X,DEG,OPTIONS) declares Pop
% on the output spaces of the 'copvar' or 'cdopvar' object X, with the
% domain of X, so that X'*Pop*X is defined. For the input spaces of X,
% pass X'.
%
% INPUT
% - prog:        'struct', the LPI program;
% - dims:        double, M x 1 positive integers, dims(k) = m_k the number
%                of components of space k; a scalar is used for every
%                space; empty for no space (POP is then 0 x 0, with no
%                decision variable, and PROG is returned unchanged);
% - spaces:      1 x M cell of 'cellstr', spaces{k} the names of the
%                variables of s^k ({} for R^{m_k}); a plain 'cellstr' is
%                one space, so several spaces are nested, as {{},{'s'}};
% - dom:         (optional) double, nv x 2, dom(d,:) the interval of
%                variable d of the registry (the sorted union of all
%                spaces), or one 1 x 2 row for every variable; or a
%                'struct' with fields vars ('cellstr') and dom (double, one
%                row per name). Omitted, as POSCOPVAR(PROG,DIMS,SPACES,DEG,
%                OPTIONS): the domain of PROG, by name;
% - X:           'copvar' or 'cdopvar' object, in place of DIMS, SPACES
%                and DOM (second form above);
% - deg:         (optional; omitted or [], the storage row of the program
%                settings: DEG OMITTED below)
%                the degrees of Z^alpha(theta,s): a scalar (int and mult),
%                or a 'struct' with the optional fields
%                  int     scalar or 1 x nv, the degree in each theta_d
%                          (default 1);
%                  mult    scalar or 1 x nv, the degree in each s_d of
%                          s^k, 0 where alpha_d = 1 (default int);
%                  joint   scalar, a cap on the total degree (default
%                          none);
%                  subset  2^(2nv) caps, one per subset of [theta, s^k],
%                          for a space with every registry variable only;
%                or a 1 x M cell of these, one per space, each of which
%                may be a cell with one per basis operator of the space;
% - options:     (optional) 'struct' with the optional fields
%                  psatz    double, the weight g: 0 (default, g = 1), 1
%                           (g = prod_d (theta_d-a_d)*(b_d-theta_d)), or
%                           2d+1 / 2d+2, d = 1..nv (g = (theta_d-a_d)/L_d
%                           or (b_d-theta_d)/L_d, L_d = b_d-a_d); any
%                           other value is an error;
%                  psatz_rows, psatz_rows_offset
%                           R x nv rows of per-direction factor codes
%                           (0 none, 1 (theta_d-a_d)/L_d, 2
%                           (b_d-theta_d)/L_d, 3 their product), one
%                           more Gram matrix per row with weight g the
%                           product of the factors of the row (face
%                           products), and their degree offsets;
%                  sep      logical, scalar or 1 x nv: where true, the
%                           lower and upper basis operators of that
%                           direction are replaced by the full integral
%                           (alpha_d = 4);
%                  include  1 x M cell, entry k selecting basis operators
%                           of space k: multi-index rows over the
%                           variables of s^k, a logical mask, or linear
%                           indices (default all).
%
% DEG OMITTED (or []): the storage row lpi.pos of the program settings
% gives DEG and the include and sep it implies, which OPTIONS may then not
% set.
% Fields other than len_R are 1 x 3 cells, entry e for a space in e
% variables, the last entry for every higher e; lengths are lengths of
% monomial vectors, 0 = absent (a pattern with a zero length is dropped):
%   len_R      scalar, the multiplier length on R^q (default 1: identity);
%   len_inc    1 x e, d1: incomplete moments, pattern c (lower integral);
%   len_com    1 x e, d2: complete moments, pattern b;
%   len_sos    3 x e, d3: multiplier lengths of the blocks a, b, c;
%   patterns   'uniform' (a..a, b..b), 'separable' ({a,b}^e), 'all'
%              ({a,b,c}^e), a cellstr of words ('aab'), or rows (1 a,
%              2 b, 3 c);
%   cap        (optional) 1 x e, Inf = none: theta_k+s_k <= cap_k in an
%              owned direction k; in a direction k the space does not
%              own, theta_k <= the first cap of the largest owner;
%   len_cross  (optional) scalar, the monomial length in a direction the
%              space does not own; absent: the largest int of the basis
%              operator in its owned directions;
%   terms, offsets  'none', []: the storage carries no term (an error
%              otherwise).
% Codes: a -> 1; c -> 2; b -> 3 in a direction where some space carries
% c, else 4 (sep true there). In an owned direction, mult = 0, d2-1,
% d1-1 and int = d3(a,:)-1, d3(b,:)-1, d3(c,:)-1 for a, b, c. Without
% len_com the fields mult, int, joint, sep of lpi.pos are read: per
% space, entry min(e,3): the multiplier basis of degree mult, integral
% bases of degree int and total joint; a direction owned by a space in
% two or more variables, or any with sep set, takes the full integral.
%
% OUTPUT
% - prog:        'struct', the program with Q declared positive
%                semidefinite and Pop appended to prog.sop_pos;
% - Pop:         'cdopvar' object, M x M blocks on X, the operator above;
% - Qcell:       N_b x N_b cell, Qcell{b,b'} the 'cellstr' names of the
%                decision variables of block (b,b') of Q, N_b the number
%                of basis operators of all spaces;
% - basis_list:  double, N_b x (1+nv), row b the basis operator b: its
%                space k, then its multi-index entry alpha_d in each
%                registry direction, 0 where s^k lacks variable d; rows
%                ordered by space, as Qcell.
%
% See also POSCOPVAR_DIRECT, COPQUADVAR, POSSOPVAR, SOPQUADVAR,
%          LPI_EQ_CDOPVAR, LPI_INEQ_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poscopvar
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
% Initial coding MMP, 09/21/2026

% Check the inputs
% DOM omitted: the 4th input is DEG (a scalar, a 'struct' without field
% dom, or a 'cell', never n x 2), and the domain is that of PROG by name,
% per element, since the varname of a polynomial vector is sorted while
% the rows of prog.dom follow the program order
nin = nargin;
% Operator form POSCOPVAR(PROG,X,DEG,OPTIONS): the output spaces of X and
% its domain stand for DIMS, SPACES, DOM; DEG and OPTIONS shift by two
if isa(dims,'copvar') || isa(dims,'cdopvar')
    narginchk(2,4);
    if nin>=4,  options = dom;  end
    if nin>=3,  deg = spaces;   end
    X = dims;   [spaces,dims] = copvar_space_list(X,'out');
    dom = struct('vars',{X.vars},'dom',X.dom);  nin = nin+2;
end
% DOM omitted, also with DEG omitted, POSCOPVAR(PROG,DIMS,SPACES): shift
% DEG and OPTIONS by one and take the domain of PROG, variable by variable
if nin==3 || (nin>=4 && ~((isnumeric(dom) && (isempty(dom) || size(dom,2)==2)) ...
               || (isstruct(dom) && isfield(dom,'dom'))))
    if nin>=5,  options = deg;  end
    if nin>=4,  deg = dom;      end
    nin = nin+1;
    v = cell(1,size(prog.dom,1));
    for k = 1:numel(v),     v{k} = prog.vartable(k).varname{1};     end
    dom = struct('vars',{v},'dom',prog.dom);
end
% DEG omitted or empty: the degrees of the program settings (below)
usedef = nin<5 || isempty(deg);
% OPTIONS: an empty struct when omitted, else a 'struct'
if nin<6 || isempty(options)
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end

% No space (an empty channel of a PIE, as the W of an H2 LPI on R^0): the
% 0 x 0 operator, with no block and no decision variable; PROG unchanged
if isempty(dims)
    meta0 = struct('vars',{cell(1,0)},'dom',zeros(0,2),'space_out',false(0,0),...
                   'space_in',false(0,0),'dim_out',zeros(0,1),'dim_in',zeros(0,1),...
                   'Zd',{cell(0,1)});
    Pop = cdopvar(cell(0,0),meta0);     Qcell = {};     basis_list = zeros(0,1);
    return
end

meta = parse_copvar_spaces(dims,spaces,dom);    % the registry, for nv
% DEG omitted: the degrees of prog.settings.lpi.pos (the struct ps below
% for a program without settings) space by space, entry min(d,3) for a space
% in d variables: on R^n the identity basis; otherwise the multiplier
% basis of degree mult and integral bases of degree int in each variable
% and joint in total (int < 0: none), the lower and upper integrals in one
% variable, one full integral in a separable direction. SEP is per
% direction: a direction that a space in two or more variables owns is
% separable for every space, any other when lpi.pos.sep is set.
if usedef
    if any(isfield(options,{'include','sep'}))
        error("With DEG omitted the settings give include and sep; OPTIONS may not set them.")
    end
    ps = struct('mult',[1 0 0],'int',[1 1 1],'joint',[2 2 2],'sep',0);
    if isfield(prog,'settings') && isfield(prog.settings,'lpi') && isfield(prog.settings.lpi,'pos')
        ps = prog.settings.lpi.pos;
    end
    % lpi.pos with len_com: the storage row (header, DEG OMITTED); else the
    % fields mult, int, joint, sep by the rule above
    if isfield(ps,'len_com')
        [deg,sepd,incl] = row_storage(ps,meta.space_out);
    else
    own = meta.space_out;   nd = sum(own,2);    e = @(x,c) x(min(c,end));
    sepd = logical(ps.sep(1)) | any(own(nd>=2,:),1);
    deg = cell(1,numel(nd));    incl = cell(1,numel(nd));
    for k = 1:numel(nd)
        % R^n: the identity basis; else the multiplier basis and the lower
        % and upper integral bases, or the full integral where separable
        if nd(k)==0,    deg{k} = struct('int',0);   continue;   end
        c = min(nd(k),3);   m = e(ps.mult,c);   b = e(ps.int,c);
        Zi = struct('int',b,'mult',b,'joint',e(ps.joint,c));
        deg{k} = {struct('int',m,'mult',0,'joint',m), Zi, Zi};
        incl{k} = [1;2;3];
        if all(sepd(own(k,:)))
            deg{k} = deg{k}(1:2);   incl{k} = [ones(1,nd(k)); 4*ones(1,nd(k))];
        end
        if b<0,     deg{k} = deg{k}(1);     incl{k} = incl{k}(1,:);     end
    end
    end
    options.sep = sepd;     options.include = incl;
end
% Check options.psatz: one code, 0, 1 or a face (theta_d-a_d) or
% (b_d-theta_d), code 2d+1 / 2d+2; a row of codes and the box quadratic
% codes 2nv+2+d are rejected
nv = numel(meta.vars);      ps = 0;
if isfield(options,'psatz') && ~isempty(options.psatz)
    ps = options.psatz;
end
if ~isscalar(ps) || ~ismember(ps,[0,1,3:2*nv+2])
    error("'psatz' should be 0, 1, or a face code 2d+1 / 2d+2 with d "...
          +"in 1..nv, here nv = "+num2str(nv)+"; 2 is "...
          +"a ball weight, not supported (1 is the box product).")
end
% Declare Q and build Pop
[prog,Pop,Qcell,info] = poscopvar_direct(prog,dims,spaces,dom,deg,options);
basis_list = info.basis_list;                   % space, then multi-index
% Append Pop to prog.sop_pos, the record of the positive operators of
% PROG; MATLAB copies on write, so the record shares the arrays of Pop
if ~isfield(prog,'sop_pos'),    prog.sop_pos = {};  end
prog.sop_pos{end+1} = Pop;

end


function [deg,sepv,incl] = row_storage(row,own)
% [DEG,SEPV,INCL] = ROW_STORAGE(ROW,OWN) the degree specifications, the
% separable directions and the include rows of the basis operators given
% by the storage row ROW (header, DEG OMITTED) on the spaces of OWN.
% INPUT  row:   'struct', lpi.pos with field len_com;
%        own:   M x nv logical, row k the registry variables of space k.
% OUTPUT deg:   1 x M cell, per space a cell of one degree struct per basis
%               operator (int, mult, cap: 1 x nv), or a struct on R^q;
%        sepv:  1 x nv logical, true where no owning space carries c;
%        incl:  1 x M cell, the multi-index rows of each space.
% Small arrays only.
[M,nv] = size(own);     nd = sum(own,2);
e_of = @(k) min(max(nd(k),1),3);    % the entry: the last for e > 3
if ~isfield(row,'len_R') || isempty(row.len_R),     row.len_R = 1;  end
% The storage carries no term
if isfield(row,'terms') && ~isempty(row.terms)
    t = row.terms;  if ~iscell(t),  t = {t};    end
    if ~all(cellfun(@(x) isempty(x) || isequal(char(x),'none'),t))
        error("The storage carries no term: lpi.pos.terms should be 'none'.")
    end
end
% Pattern rows per space (1 a, 2 b, 3 c); a zero length drops the pattern
pats = cell(1,M);
for k = 1:M
    if nd(k)==0,    continue;   end
    e = e_of(k);    n = nd(k);
    % Pattern set: a named set, words over a, b, c, or rows of codes
    pt = row.patterns{e};
    if ischar(pt)
        switch pt
            case 'uniform',     A = [ones(1,n); 2*ones(1,n)];
            case 'separable',   A = multiindex_grid(repmat({[1 2]},1,n));
            case 'all',         A = multiindex_grid(repmat({[1 2 3]},1,n));
            otherwise,  error("Unknown pattern set '"+pt+"' in lpi.pos.patterns.")
        end
    elseif iscell(pt)   % words over a, b, c
        A = zeros(numel(pt),n);
        for i = 1:numel(pt),    w = pt{i};    A(i,:) = (w=='b')*1 + (w=='c')*2 + 1;   end
    else
        A = pt;
    end
    % Lengths of this entry; drop a pattern with a zero length in any
    % direction
    d1 = fit_len(row.len_inc{e},n);   d2 = fit_len(row.len_com{e},n);   d3 = row.len_sos{e};
    if size(d3,2)~=n,   d3 = repmat(d3(:,1),1,n);    end
    keep = true(size(A,1),1);
    for i = 1:size(A,1)
        for d = 1:n
            a = A(i,d);
            if (a==2 && d2(d)==0) || (a==3 && d1(d)==0) || d3(a,d)==0,   keep(i) = false;  end
        end
    end
    pats{k} = A(keep,:);
end
% Separable directions: no owning space carries c there (b -> 4, the full
% integral, spans b and c with equal multiplier lengths)
sepv = true(1,nv);
for k = 1:M
    if nd(k)==0,    continue;   end
    idx = find(own(k,:));
    for d = 1:nd(k)
        if any(pats{k}(:,d)==3),    sepv(idx(d)) = false;    end
    end
end
% The number of variables of the largest space owning each direction
big = zeros(1,nv);
for d = 1:nv,   o = find(own(:,d));  if ~isempty(o),  big(d) = max(nd(o));  end,  end
% Degrees and include rows
deg = cell(1,M);    incl = cell(1,M);
for k = 1:M
    % R^q: monomials of length len_R; else one basis operator per pattern
    if nd(k)==0,    deg{k} = struct('int',row.len_R-1);    continue;   end
    e = e_of(k);    n = nd(k);  idx = find(own(k,:));
    d1 = fit_len(row.len_inc{e},n);   d2 = fit_len(row.len_com{e},n);   d3 = row.len_sos{e};
    if size(d3,2)~=n,   d3 = repmat(d3(:,1),1,n);    end
    A = pats{k};    rows = zeros(size(A));   dk = cell(1,size(A,1));
    for i = 1:size(A,1)
        % Owned directions: the code and the degrees in s (mult) and
        % theta (int) of the pattern letter
        mu = zeros(1,nv);   in = zeros(1,nv);
        for d = 1:n
            a = A(i,d);     g = idx(d);
            switch a    % b: 4 (full) where separable, else 3 (upper)
                case 1,     rows(i,d) = 1;                      mu(g) = 0;
                case 2,     rows(i,d) = 4 - ~sepv(g);           mu(g) = d2(d)-1;
                case 3,     rows(i,d) = 2;                      mu(g) = d1(d)-1;
            end
            in(g) = d3(a,d)-1;
        end
        % A direction the space does not own: len_cross - 1 when given,
        % else the largest owned int of the basis operator
        for g = find(~own(k,:))
            in(g) = max(in(idx));
            if isfield(row,'len_cross') && ~isempty(row.len_cross) && ~isempty(row.len_cross{e})
                in(g) = row.len_cross{e} - 1;
            end
        end
        dk{i} = struct('int',in,'mult',mu);
        % Per-direction cap: owned directions from this entry,
        % a non-owned one from the entry of its largest owner
        if isfield(row,'cap') && ~isempty(row.cap)
            cv = Inf(1,nv);
            ce = fit_len(row.cap{e},n);
            for d = 1:n,    cv(idx(d)) = ce(d);     end
            for g = find(~own(k,:))
                cb = row.cap{min(big(g),3)};    cv(g) = cb(1);
            end
            dk{i}.cap = cv;
        end
    end
    deg{k} = dk;    incl{k} = rows;
end
end


function v = fit_len(v,n)
% V = FIT_LEN(V,N) a row of N lengths from a scalar or a row of N.
% INPUT  v:  double, a scalar or N entries;   n: double, the row length.
% OUTPUT v:  double, 1 x N; an error for any other number of entries.
v = reshape(v,1,[]);
if isscalar(v),     v = repmat(v,1,n);
elseif numel(v)~=n
    error("A length row of lpi.pos has "+numel(v)+" entries for a space in "+n+" variables.")
end
end
