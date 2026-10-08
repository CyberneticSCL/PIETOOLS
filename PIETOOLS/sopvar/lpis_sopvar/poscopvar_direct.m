function [prog,Pop,Qcell,info] = poscopvar_direct(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,INFO] = POSCOPVAR_DIRECT(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares the positive semidefinite, self-adjoint 'cdopvar' decision
% operator of 'poscopvar' on a concatenation of mixed spaces by the DIRECT
% MAP from the Gram matrix to the kernel coefficients: no polynomial
% product, no coefficient sheet and no operator composition in between.
%
% THE MAP. With the basis operators of 'copquadvar' (Sec. 9.1 of the sopvar
% document), Gram row u = (component c, matrix row r, monomial theta^p s^j)
% and column v = (c', t, theta^q s'^l), and one Positivstellensatz weight
% g(theta) = prod_d g_d(theta_d), the kernel of Pop is
%
%   sum_{u,v} Q_uv s^j s'^l e_r e_t' prod_d int I_{alpha_d}(theta_d - s_d)
%                      I_{beta_d}(theta_d - s'_d) g_d(theta_d) theta_d^{p_d+q_d} dtheta_d,
%
% alpha, beta the multi-indices of c and c' (1 multiplier, 2 lower, 3 upper,
% 4 full). Every factor is an entry of the 1-D table of Sec. 6.1 with
% G = g_d theta^n, n = p_d + q_d, and calG its antiderivative (lower cell
% s' <= s, upper cell s <= s'):
%
%   (1,1) mult: G(s)        (1,2) low: G(s)            (1,3) up: G(s)
%   (2,1) up:   G(s')       (3,1) low: G(s')
%   (2,2) low:  calG(b)-calG(s)      up: calG(b)-calG(s')
%   (3,3) low:  calG(s')-calG(a)     up: calG(s)-calG(a)
%   (2,3) up:   calG(s')-calG(s)     (3,2) low: calG(s)-calG(s')
%   a full index is the sum of the two one-sided rows (or columns),
%
% and a direction only one space has, or neither has, is the one-sided or
% definite integral of the same G. Each entry is a short list of (left
% exponent, right exponent, coefficient); the N-D entry is the Kronecker
% product of the per-direction lists, the lift monomials add j and l to the
% exponents, and every term becomes a triplet (Gram position, coefficient
% column of the cell, value). One 'sparse' per cell, then the symmetric
% placement of the decision variables, gives params.B of every block. The
% bases ZL, ZR are the exponents reached, so no padding is stored.
%
% TWO PATHS (OPTIONS.path):
%   'general'  the tables are built for the term's own weight g_d, given by
%              its coefficients, at the Gram's own degrees;
%   'fast'     (default) a canonical weight is a product of degree-1
%              factors per direction, so it moves onto the Gram side as a
%              shift: g Z_w'^T S Z_w' = Z_w^T (E^T S F) Z_w with
%              (theta-a) v_{w'} = E_a v_w and (b-theta) v_{w'} = E_b v_w, the
%              product term on both sides, a face term on the left only.
%              The tables are then those of g = 1 (single monomials) on the
%              shifted index set, and the position map, which depends on the
%              degrees, domain and spaces but not on the program, is CACHED
%              across calls (OPTIONS.cache) and reused by every term that
%              lands on the same index set, e.g. the plain term at w and the
%              product term at w-1 of 'lpi_ineq'. The shift route is taken
%              only when its map is cached (INFO.route 'fast-shift'); on a
%              miss the term's own map is built and cached ('fast-own'),
%              since the shifted map is the plain map on the enlarged index
%              set with the shift's fan-out on top.
% Both paths give the same Pop up to rounding; measured in
% test_poscopvar_direct.
%
% INPUT
% - prog, dims, spaces, dom: as 'poscopvar';
% - deg:     as 'copquadvar' (int, mult, joint, subset; a scalar; a cell per
%            space, each a cell per basis operator). The aliases 'lift'
%            (for mult) and 'weight' (for int) of 'poscopvar_lift' are
%            accepted. Joint and subset caps ARE supported: the map needs
%            only the index set of the Gram;
% - options: (optional) struct with fields
%   sep, include   as 'copquadvar';
%   psatz          term codes as 'poscopvar', one or a ROW of them: 0
%                  (default), 1 (the box product), 2d+1 / 2d+2 (lower /
%                  upper face in direction d, normalised by the length).
%                  With several codes Pop is the SUM of the terms, one Gram
%                  each, declared in the given order, all on one decision
%                  list and built in one pass: the multipliers are added
%                  before the lift is applied;
%   psatz_offset   per term (or one for all, default 0): the term's 'int'
%                  lowered by it in every direction, floored at 0, 'mult'
%                  kept, a 'joint' cap lowered by it per direction the
%                  weight acts in ('subset' caps refused). [0 1] at codes
%                  [0 1] is the Markov-Lukacs pair S0 + g S1 of 1-D, which
%                  the plain term's cached map serves at once;
%   path           'fast' (default) or 'general', see above;
%   cache          logical, default true: keep the position maps of this
%                  session. poscopvar_direct('clear') empties the cache;
%
% OUTPUT
% - prog:    the program with the Gram declared ('sosquadvar', type 'pos');
% - Pop:     M x M 'cdopvar' on the given spaces, Pop = Pop' >= 0, on the
%            same decision variables, in the same order, as 'poscopvar'
%            would declare for the same inputs on the same program;
% - Qcell:   N_b x N_b cell, Qcell{c,c'} the names of the decision
%            variables of the Gram block of basis operators c and c', as
%            'copquadvar' returns them; with several terms a 1 x nt cell of
%            such cells, one per term;
% - info:    struct with fields basis_list, E (the monomial grid of every
%            basis operator, first term), N (Gram dimensions, one per
%            term), ndec (all terms), path, code, offset, route
%            ('fast-shift', 'fast-own' or 'general'; a cell with several
%            terms), cache_hit, Qdpvar (the Gram as the N x N 'dpvar' of
%            'sosquadvar', for 'sosgetsol'; a cell with several terms),
%            terms (per-term struct array: code, offset, E, N, off, Nc,
%            names, ndec, route, cache_hit, times) and the summed times
%            t_declare, t_map, t_product, t_blocks, t_total in seconds.
%
% NOTES
% Cost: the Gram declaration is 'sosquadvar' with the constant basis, O(N^2)
% in the Gram dimension N; the map is one pass over the Gram positions with
% a fan-out equal to the product of the per-direction list lengths, which is
% the nnz of B plus the duplicates 'sparse' sums. Nothing is proportional
% to the number of decision variables except B itself and the symmetric
% placement D, which has two entries per variable.
%
% See also POSCOPVAR, COPQUADVAR, POSCOPVAR_LIFT, LIFT_COPVAR,
% PRUNE_ZERO_MONOMIALS.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poscopvar_direct
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
% Initial coding MMP, 10/07/2026. The direct Gram-to-kernel map proposed
%                after profiling 'copquadvar' (57% of its time in the
%                symbolic product, the sheet conversion and the per-pair
%                term objects) and the separated form (three operator
%                objects on union bases, 34 GB in 3-D).
% MMP, 10/08/2026: Several Positivstellensatz terms in one call ('psatz' a
%                row, 'psatz_offset' per term): one Gram per term, the maps
%                shared where the index sets coincide, every term's map
%                embedded in the union of the bases and summed on one
%                decision list, so a sum of terms costs one block build and
%                no container addition. For 'lpi_ineq_sop'. One term is
%                unchanged (the same program, Qcell and info as before).

persistent CACHE
if nargin>=1 && (ischar(prog) || isstring(prog))
    if strcmpi(char(prog),'clear')
        CACHE = [];     prog = [];  Pop = [];  Qcell = {};  info = struct('cleared',true);
        return
    end
    error("The first input is the program, or 'clear' to empty the cache.")
end
if nargin<5
    error("Not enough input arguments.")
end
if nargin<6 || isempty(options),    options = struct(); end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if isfield(options,'type') && ~isempty(options.type) && ~strcmp(char(options.type),'pos')
    error("'poscopvar_direct' declares a positive operator (type 'pos').")
end
if isempty(CACHE),  CACHE = containers.Map('KeyType','char','ValueType','any');    end
codes = 0;
if isfield(options,'psatz') && ~isempty(options.psatz),  codes = double(options.psatz(:).');   end
offs = zeros(size(codes));
if isfield(options,'psatz_offset') && ~isempty(options.psatz_offset)
    offs = double(options.psatz_offset(:).');
    if isscalar(offs),  offs = repmat(offs,size(codes));
    elseif numel(offs)~=numel(codes)
        error("'psatz_offset' should be a scalar or have one entry per psatz term.")
    end
end
if any(offs<0 | offs~=round(offs))
    error("'psatz_offset' should be nonnegative integers.")
end
nt = numel(codes);
path = 'fast';
if isfield(options,'path') && ~isempty(options.path),    path = lower(char(options.path));   end
if ~any(strcmp(path,{'fast','general'}))
    error("'path' should be 'fast' or 'general'.")
end
usecache = true;
if isfield(options,'cache') && ~isempty(options.cache),  usecache = logical(options.cache);  end
t_all = tic;

% % % Spaces, registry, components: 'copquadvar's enumeration, through
% 'lift_copvar' (multiindex_grid, direction 1 fastest, then include and sep).
[meta,sp] = parse_copvar_spaces(dims,spaces,dom);
M = numel(sp);  vars = meta.vars;   nv = numel(vars);   dom = meta.dom;
mask = meta.space_out;  mk = meta.dim_out(:);
own = cell(1,M);
for k = 1:M,    own{k} = reshape(find(mask(k,:)),1,[]);     end
if any(~ismember(codes,[0,1,3:2*nv+2]))
    error("'psatz' entries should be 0, 1, or face codes 2d+1 / 2d+2 with d in 1..nv.")
end
lopt = struct();
if isfield(options,'sep'),      lopt.sep = options.sep;         end
if isfield(options,'include'),  lopt.include = options.include; end
[~,li] = lift_copvar(dims,spaces,dom,1,lopt);
basis_list = li.basis_list;     Nb = size(basis_list,1);    bsp = basis_list(:,1);
bix = zeros(Nb,1);  nbk = zeros(1,M);
for c = 1:Nb,   nbk(bsp(c)) = nbk(bsp(c))+1;    bix(c) = nbk(bsp(c));   end

% % % The degree specification per space (or per basis operator).
deg = alias_deg(deg);
if iscell(deg)
    if numel(deg)~=M,   error("A cell 'deg' should have one entry per space.");  end
    deg_sp = reshape(deg,1,[]);
else
    deg_sp = repmat({deg},1,M);
end

% % % Per term: the Gram's monomial grids ('int' lowered by the term's
% offset), its declaration, its weight and position map. The maps are
% applied once the common decision list of all terms is known.
terms = struct('code',{},'offset',{},'E',{},'N',{},'off',{},'Nc',{},'names',{},'ndec',{},...
               'route',{},'cache_hit',{},'Qdpvar',{},'PM',{},'SE',{},'SF',{},'useshift',{},...
               't_declare',{},'t_map',{},'t_product',{});
for c = 1:nt
    % The monomial grid of every basis operator, as copquadvar's
    % monomial_bases: exponents over [theta_1..theta_nv, s^k], 'mult' zero
    % in a multiplier direction (delta identifies s_d with theta_d).
    E = cell(1,Nb);
    for b = 1:Nb
        k = bsp(b);     spec = deg_sp{k};
        if iscell(spec)
            if numel(spec)~=nbk(k)
                error("A cell degree specification for space "+num2str(k)...
                      +" should have one entry per basis operator of that space.")
            end
            spec = spec{bix(b)};
        end
        spec = reduce_spec(spec,offs(c),ng_of(codes(c),nv),nv);
        d = process_degrees_one(spec,nv);
        ok = own{k};
        if ~isempty(d.subset) && numel(ok)~=nv
            error("A 'subset' degree array is only accepted for a space holding "...
                  +"every registry variable; space "+num2str(k)+" does not.")
        end
        caps_mult = d.mult(ok);
        caps_mult(basis_list(b,1+ok)==1) = 0;
        E{b} = build_exponent_grid([d.int,caps_mult],d.joint,d.subset);
    end

    % The Gram: rows (basis operator outer, matrix row, monomial inner), the
    % index of 'sosquadvar', declared with the constant basis so that the
    % names and the lower-triangular numbering are those 'copquadvar' gets.
    Tc = cellfun(@(e) size(e,1),E);     Tc = Tc(:);
    Nc = Tc.*mk(bsp);
    off = [0; cumsum(Nc)];  N = off(end);
    t0 = tic;
    [prog,Pq] = sosquadvar(prog,{polynomial(1)},{polynomial(1)},N,N,'pos');
    names = Pq{1,1}.dvarname;
    if ~iscellstr(names),   names = cellstr(string(names));    end
    names = reshape(names,[],1);
    ndec_c = numel(names);
    if ndec_c~=N*(N+1)/2
        error("Internal error: the Gram declaration returned %d variables for N = %d.",ndec_c,N)
    end
    t_declare = toc(t0);

    % The term's weight: general path, the coefficients of g_d; fast path,
    % g = 1 with the degree-1 factors as shifts of the Gram index.
    t0 = tic;
    useshift = strcmp(path,'fast') && codes(c)~=0;
    [gco,shL,shR] = term_weight(codes(c),nv,dom,path);
    [ELt,ERt,SE,SF] = target_grids(E,nv,shL,shR);
    key = cache_key(basis_list,ELt,ERt,mk,dom,gco,own);
    hit = usecache && isKey(CACHE,key);
    if useshift && ~hit
        % A cold shifted map is the plain map on the ENLARGED index set,
        % with the shift's fan-out on top (measured 3-D, product at int 1:
        % 26 s against 17 s for the term's own map), so the shift route is
        % taken only when its map is already cached; otherwise the term's
        % own map is built and cached under its own key.
        useshift = false;
        [gco,shL,shR] = term_weight(codes(c),nv,dom,'general');
        [ELt,ERt,SE,SF] = target_grids(E,nv,shL,shR);
        key = cache_key(basis_list,ELt,ERt,mk,dom,gco,own);
        hit = usecache && isKey(CACHE,key);
    end
    if useshift,                        route = 'fast-shift';
    elseif strcmp(path,'fast'),         route = 'fast-own';
    else,                               route = 'general';
    end
    % The position map: Gram positions (u,v) of the target index set to the
    % coefficient columns of every block and cell. Cached on everything it
    % depends on.
    if hit
        PM = CACHE(key);
    else
        PM = pos_map(ELt,ERt,basis_list,mk,own,dom,gco,nv);
        if usecache,    CACHE(key) = PM;    end
    end
    t_map = toc(t0);
    terms(c) = struct('code',codes(c),'offset',offs(c),'E',{E},'N',N,'off',off,'Nc',Nc,...
                      'names',{names},'ndec',ndec_c,'route',route,'cache_hit',hit,...
                      'Qdpvar',Pq{1,1},'PM',PM,'SE',{SE},'SF',{SF},'useshift',useshift,...
                      't_declare',t_declare,'t_map',t_map,'t_product',0);
end

% % % The decision variables of all terms in one list, sorted as copquadvar
% sorts its names, and each term's placement: dvar k sits at (u,v) and
% (v,u) of its Gram, rows in the common list; the shifts move the source
% positions to the target positions.
allnames = vertcat(terms.names);
[Zd,~,rowall] = unique(allnames);
ndec = numel(Zd);
D = cell(1,nt);     pos = 0;
for c = 1:nt
    t0 = tic;
    rm = rowall(pos+(1:terms(c).ndec));     pos = pos+terms(c).ndec;
    N = terms(c).N;
    [U,V] = ndgrid(1:N,1:N);
    kk = lind_sym(U(:),V(:),N);
    D{c} = sparse(rm(kk),U(:)+(V(:)-1)*N,1,ndec,N^2);
    if terms(c).useshift
        SEf = blkdiag_comp(terms(c).SE,mk(bsp));    SFf = blkdiag_comp(terms(c).SF,mk(bsp));
        D{c} = D{c}*kron(SFf,SEf);
    end
    terms(c).t_product = toc(t0);
end

% % % The blocks: every term's map, embedded in the union of the bases the
% terms reach, summed over the terms (disjoint rows).
t0 = tic;
Cblk = cell(M,M);
for k = 1:M
    for kp = 1:M
        ok = own{k};    op = own{kp};
        ZL = cell(1,numel(ok));     ZR = cell(1,numel(op));
        for c = 1:nt
            b = terms(c).PM.blk{k,kp};
            for i = 1:numel(ok),    ZL{i} = union(ZL{i},b.ZL{i});   ZL{i} = ZL{i}(:);  end
            for i = 1:numel(op),    ZR{i} = union(ZR{i},b.ZR{i});   ZR{i} = ZR{i}(:);  end
        end
        NL = prod([cellfun(@numel,ZL),1]);  NR = prod([cellfun(@numel,ZR),1]);
        nAB = mk(k)*NL*mk(kp)*NR;
        n3 = terms(1).PM.blk{k,kp}.n3;      ncell = 3^n3;
        Bsum = cell(ncell,1);
        for c = 1:nt
            b = terms(c).PM.blk{k,kp};
            Ec = embed_cols(b.ZL,b.ZR,ZL,ZR,mk(k),mk(kp));
            for q = 1:ncell
                if isempty(b.Bpos{q}),  continue,   end
                Bq = D{c}*b.Bpos{q};
                if ~isempty(Ec),    Bq = Bq*Ec;     end
                if isempty(Bsum{q}),    Bsum{q} = Bq;   else,   Bsum{q} = Bsum{q}+Bq;  end
            end
        end
        params = struct('A',{cell(ncell,1)},'B',{cell(ncell,1)});
        for q = 1:ncell
            params.A{q} = sparse(nAB,1);
            if isempty(Bsum{q}),    params.B{q} = sparse(ndec,nAB);
            else,                   params.B{q} = Bsum{q};
            end
        end
        params.A = reshape(params.A,[3*ones(1,n3),1,1]);
        params.B = reshape(params.B,[3*ones(1,n3),1,1]);
        vio = struct('out',{vars(ok)},'in',{vars(op)});
        dio = struct('out',dom(ok,:),'in',dom(op,:));
        Cblk{k,kp} = sdopvar(params,vio,Zd,ZL,ZR,dio,[mk(k),mk(kp)]);
    end
end
Pop = cdopvar(Cblk);
t_blocks = toc(t0);

% % % Qcell: the names of the Gram block of every pair of basis operators,
% per term (one cell when there is one term).
Qcell = {};
if nargout>=3
    Qc = cell(1,nt);
    for c = 1:nt
        off = terms(c).off;     Nc = terms(c).Nc;   N = terms(c).N;    names = terms(c).names;
        Qc{c} = cell(Nb,Nb);
        for b = 1:Nb
            for bp = 1:Nb
                [Uc,Vc] = ndgrid(off(b)+(1:Nc(b)),off(bp)+(1:Nc(bp)));
                Qc{c}{b,bp} = reshape(names(lind_sym(Uc(:),Vc(:),N)),size(Uc));
            end
        end
    end
    if nt==1,   Qcell = Qc{1};  else,   Qcell = Qc;     end
end
if nt==1,   route = terms(1).route;     Qd = terms(1).Qdpvar;
else,       route = {terms.route};      Qd = {terms.Qdpvar};
end
info = struct('basis_list',basis_list,'E',{terms(1).E},'N',[terms.N],'ndec',ndec,'path',path,...
              'route',{route},'code',codes,'offset',offs,'cache_hit',[terms.cache_hit],...
              't_declare',sum([terms.t_declare]),'t_map',sum([terms.t_map]),...
              't_product',sum([terms.t_product]),'t_blocks',t_blocks,'t_total',toc(t_all),...
              'Qdpvar',{Qd},'terms',rmfield(terms,{'PM','SE','SF','useshift','Qdpvar'}));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = ng_of(code,nv)
% Directions in which the weight of a term has positive degree.
if code==0,         n = 0;
elseif code==1,     n = nv;
else,               n = 1;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = reduce_spec(spec,off,ng,nv)
% The degree specification of a Psatz term: 'int' lowered by OFF in every
% direction, floored at 0, 'mult' kept at what it was (its default is the
% ORIGINAL 'int'), a 'joint' cap lowered by OFF per direction the weight
% acts in. A 'subset' cap cannot move with 'int' and is refused.
if off==0,  return,     end
if isnumeric(spec) && isscalar(spec)
    spec = struct('int',spec,'mult',spec);
elseif ~isstruct(spec)
    error("Monomial degrees should be specified as a scalar or 'struct' object.")
end
if isfield(spec,'subset') && ~isempty(spec.subset)
    error("A 'psatz_offset' is not defined for a 'subset' degree specification.")
end
w = 1;
if isfield(spec,'int') && ~isempty(spec.int),   w = spec.int;   end
if ~isfield(spec,'mult') || isempty(spec.mult), spec.mult = w;  end
w = reshape(w,1,[]);    if isscalar(w),   w = repmat(w,1,max(nv,1));    end
spec.int = max(w-off,0);
if isfield(spec,'joint') && ~isempty(spec.joint)
    spec.joint = max(spec.joint-off*ng,0);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ec = embed_cols(ZLc,ZRc,ZL,ZR,m,n)
% The sparse map from a term's coefficient columns, on its own bases ZLc,
% ZRc, to the columns on the union bases ZL, ZR of the block; [] when the
% bases coincide. Columns are (row r, left monomial, column t, right
% monomial), left and right monomials first variable slowest.
if isequal(ZLc,ZL) && isequal(ZRc,ZR),  Ec = [];  return,     end
ELc = multiindex_grid(ZLc,'first_slowest');     NLc = size(ELc,1);
ERc = multiindex_grid(ZRc,'first_slowest');     NRc = size(ERc,1);
nZL = cellfun(@numel,ZL);   nZR = cellfun(@numel,ZR);
NL = prod([nZL,1]);         NR = prod([nZR,1]);
lidx = ones(NLc,1);
for i = 1:numel(ZL)
    [~,p] = ismember(ELc(:,i),ZL{i});
    lidx = lidx+(p-1)*prod([nZL(i+1:end),1]);
end
ridx = ones(NRc,1);
for i = 1:numel(ZR)
    [~,p] = ismember(ERc(:,i),ZR{i});
    ridx = ridx+(p-1)*prod([nZR(i+1:end),1]);
end
nABc = m*NLc*n*NRc;
j = (1:nABc)'-1;
rowpart = mod(j,m*NLc);     colpart = floor(j/(m*NLc));
lc = mod(rowpart,NLc)+1;    r = floor(rowpart/NLc)+1;
rc = mod(colpart,NRc)+1;    t = floor(colpart/NRc)+1;
ju = (r-1)*NL+lidx(lc)+((t-1)*NR+ridx(rc)-1)*(m*NL);
Ec = sparse(1:nABc,ju,1,nABc,m*NL*n*NR);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function k = lind_sym(u,v,N)
% The decision variable of Gram position (u,v): sosquadvar's numbering of
% the lower triangle by columns, mirrored.
U = max(u,v);   V = min(u,v);
k = U + (V-1)*N - V.*(V-1)/2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = blkdiag_comp(Sc,m)
% The shift of the whole Gram index from the per-component shifts: rows
% (component, matrix row, monomial), so component c contributes
% kron(I_m, Sc{c}).
B = cell(1,numel(Sc));
for c = 1:numel(Sc),    B{c} = kron(speye(m(c)),Sc{c});   end
S = blkdiag(B{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = alias_deg(deg)
% 'lift' and 'weight' of poscopvar_lift read as copquadvar's 'mult', 'int'.
if iscell(deg)
    for i = 1:numel(deg),   deg{i} = alias_deg(deg{i});   end
elseif isstruct(deg)
    if isfield(deg,'lift') && ~isempty(deg.lift)
        if ~isfield(deg,'mult') || isempty(deg.mult),   deg.mult = deg.lift;    end
        deg = rmfield(deg,'lift');
    end
    if isfield(deg,'weight') && ~isempty(deg.weight)
        if ~isfield(deg,'int') || isempty(deg.int),     deg.int = deg.weight;   end
        deg = rmfield(deg,'weight');
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gco,shL,shR] = term_weight(code,nv,dom,path)
% gco{d}: coefficients of g_d (constant first) for the general path (1 on
% the fast path); shL{d}, shR{d}: rows [shift, coefficient] of the left and
% right degree-1 factors on the fast path (identity on the general path).
gco = repmat({1},1,nv);     shL = repmat({[0,1]},1,nv);   shR = shL;
if code==1
    for d = 1:nv
        a = dom(d,1);   b = dom(d,2);
        if strcmp(path,'fast')
            shL{d} = [0,-a; 1,1];       shR{d} = [0,b; 1,-1];
        else
            gco{d} = [-a*b, a+b, -1];
        end
    end
elseif code>=3
    d = floor((code-1)/2);  a = dom(d,1);   b = dom(d,2);   L = b-a;
    if mod(code,2)==1,  fl = [-a,1]/L;  else,   fl = [b,-1]/L;    end
    if strcmp(path,'fast'),     shL{d} = [0,fl(1); 1,fl(2)];
    else,                       gco{d} = fl;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Et,S] = shift_grid(E,nv,sh)
% The target grid of a component under the per-direction shifts sh{d}
% (rows [e, coef]) and the sparse map S, S(t_source,t_target) = the product
% of the coefficients of the shift combination taking one to the other.
T = size(E,1);
combos = multiindex_grid(cellfun(@(s) 1:size(s,1),sh,'UniformOutput',false));
nc = size(combos,1);
Eall = zeros(nc*T,size(E,2));   rows = zeros(nc*T,1);    vals = zeros(nc*T,1);
for i = 1:nc
    ev = zeros(1,nv);   cf = 1;
    for d = 1:nv
        ev(d) = sh{d}(combos(i,d),1);   cf = cf*sh{d}(combos(i,d),2);
    end
    Es = E;     Es(:,1:nv) = Es(:,1:nv)+ev;
    Eall((i-1)*T+(1:T),:) = Es;     rows((i-1)*T+(1:T)) = (1:T)';   vals((i-1)*T+(1:T)) = cf;
end
[Et,~,ic] = unique(Eall,'rows');
S = sparse(rows,ic,vals,T,size(Et,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = cache_key(basis_list,ELt,ERt,mk,dom,gco,own)
s = jsonencode(struct('bl',basis_list,'EL',{ELt},'ER',{ERt},'m',mk,'dom',dom,...
                      'g',{gco},'own',{own}));
key = s;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = dir_list(kind,al,be,n,g,a,b)
% The 1-D table of the header for G = g(theta) theta^n: rows
% [gamma, left exponent, right exponent, coefficient]; gamma is 1 when the
% direction is not shared. kind: 3 shared, 2 output space only, 1 input
% space only, 0 neither.
K = numel(g)-1;     kk = (0:K)';     gk = reshape(g,[],1);
o = ones(K+1,1);    z = zeros(K+1,1);
Gs  = [o, n+kk, z, gk];                     % G(s)
Gsp = [o, z, n+kk, gk];                     % G(s')
As  = [o, n+kk+1, z, gk./(n+kk+1)];         % calG(s)
Asp = [o, z, n+kk+1, gk./(n+kk+1)];         % calG(s')
Cb  = [1, 0, 0, sum(gk.*b.^(n+kk+1)./(n+kk+1))];
Ca  = [1, 0, 0, sum(gk.*a.^(n+kk+1)./(n+kk+1))];
neg = @(X) [X(:,1:3), -X(:,4)];
gam = @(X,gm) [gm*ones(size(X,1),1), X(:,2:4)];
switch kind
    case 0
        L = [Cb; neg(Ca)];
    case 2
        switch al
            case 1,     L = Gs;
            case 2,     L = [Cb; neg(As)];
            case 3,     L = [As; neg(Ca)];
            case 4,     L = [Cb; neg(Ca)];
        end
    case 1
        switch be
            case 1,     L = Gsp;
            case 2,     L = [Cb; neg(Asp)];
            case 3,     L = [Asp; neg(Ca)];
            case 4,     L = [Cb; neg(Ca)];
        end
    case 3
        switch 10*al+be
            case 11,    L = gam(Gs,1);
            case 12,    L = gam(Gs,2);
            case 13,    L = gam(Gs,3);
            case 14,    L = [gam(Gs,2); gam(Gs,3)];
            case 21,    L = gam(Gsp,3);
            case 31,    L = gam(Gsp,2);
            case 41,    L = [gam(Gsp,2); gam(Gsp,3)];
            case 22,    L = [gam([Cb; neg(As)],2); gam([Cb; neg(Asp)],3)];
            case 33,    L = [gam([Asp; neg(Ca)],2); gam([As; neg(Ca)],3)];
            case 23,    L = gam([Asp; neg(As)],3);
            case 32,    L = gam([As; neg(Asp)],2);
            case 42,    L = [gam([Cb; neg(Asp)],2); gam([Cb; neg(Asp)],3)];
            case 24,    L = [gam([Cb; neg(As)],2); gam([Cb; neg(As)],3)];
            case 43,    L = [gam([Asp; neg(Ca)],2); gam([Asp; neg(Ca)],3)];
            case 34,    L = [gam([As; neg(Ca)],2); gam([As; neg(Ca)],3)];
            case 44,    L = [gam([Cb; neg(Ca)],2); gam([Cb; neg(Ca)],3)];
            otherwise,  error("Internal error: pair (%d,%d).",al,be)
        end
    otherwise
        error("Internal error: kind %d.",kind)
end
L = L(L(:,4)~=0,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PM = pos_map(ELt,ERt,basis_list,mk,own,dom,gco,nv)
% The map from the target Gram positions to the coefficient columns of
% every block and parameter cell: PM.blk{k,kp} with fields Bpos (cell over
% the 3^n3 cells, each sparse NLt*NRt x nAB), ZL, ZR, nAB, n3.
Nb = size(basis_list,1);    bsp = basis_list(:,1);     M = numel(own);
TL = cellfun(@(e) size(e,1),ELt);   TL = TL(:);
TR = cellfun(@(e) size(e,1),ERt);   TR = TR(:);
NLc = TL.*mk(bsp);      NRc = TR.*mk(bsp);
offL = [0; cumsum(NLc)];    offR = [0; cumsum(NRc)];
NLt = offL(end);    NRt = offR(end);

% Degree ranges per component and direction: theta on each side, lift.
PL = zeros(Nb,nv);  PR = zeros(Nb,nv);  JL = zeros(Nb,nv);  JR = zeros(Nb,nv);
stL = cell(1,Nb);   stR = cell(1,Nb);   tlL = cell(1,Nb);   tlR = cell(1,Nb);
for c = 1:Nb
    ok = own{bsp(c)};
    PL(c,:) = max(ELt{c}(:,1:nv),[],1);     PR(c,:) = max(ERt{c}(:,1:nv),[],1);
    if ~isempty(ok)
        JL(c,ok) = max(ELt{c}(:,nv+1:end),[],1);
        JR(c,ok) = max(ERt{c}(:,nv+1:end),[],1);
    end
    % monomial index on the target grid from the exponent row, 0 if absent
    rng = [PL(c,:)+1, JL(c,ok)+1];      stL{c} = cumprod([1,rng(1:end-1)]);
    tlL{c} = zeros(prod(rng),1);        tlL{c}(1+ELt{c}*stL{c}(:)) = (1:TL(c))';
    rng = [PR(c,:)+1, JR(c,ok)+1];      stR{c} = cumprod([1,rng(1:end-1)]);
    tlR{c} = zeros(prod(rng),1);        tlR{c}(1+ERt{c}*stR{c}(:)) = (1:TR(c))';
end
Kg = cellfun(@(g) numel(g)-1,gco);
ELmax = max(JL,[],1)+max(PL,[],1)+max(PR,[],1)+Kg+1;    % exponent ranges per direction
ERmax = max(JR,[],1)+max(PL,[],1)+max(PR,[],1)+Kg+1;
nmax = max(PL,[],1)+max(PR,[],1);

% The 1-D lists, memoized per direction on (kind, alpha, beta, n).
LST = cell(1,nv);
for d = 1:nv,   LST{d} = cell(4,5,5,nmax(d)+1);    end
    function L = getlist(d,kind,al,be,n)
        L = LST{d}{kind+1,al+1,be+1,n+1};
        if isempty(L)
            L = dir_list(kind,al,be,n,gco{d},dom(d,1),dom(d,2));
            if isempty(L),  L = zeros(0,4);     end
            LST{d}{kind+1,al+1,be+1,n+1} = L;
        end
    end

blk = cell(M,M);
for k = 1:M
  for kp = 1:M
    ok = own{k};    op = own{kp};
    D3 = intersect(ok,op);  n3 = numel(D3);   ncell = 3^n3;
    kind = zeros(1,nv);
    kind(ismember(1:nv,ok) & ismember(1:nv,op)) = 3;
    kind(ismember(1:nv,ok) & ~ismember(1:nv,op)) = 2;
    kind(~ismember(1:nv,ok) & ismember(1:nv,op)) = 1;
    cs = reshape(find(bsp==k),1,[]);    cps = reshape(find(bsp==kp),1,[]);
    % Pass 1: the per-direction matrices of every pair, and the exponents hit.
    Td = cell(numel(cs),numel(cps),nv);
    hitL = cell(1,nv);  hitR = cell(1,nv);
    for d = 1:nv,   hitL{d} = false(ELmax(d)+1,1);  hitR{d} = false(ERmax(d)+1,1);  end
    for ic = 1:numel(cs)
      c = cs(ic);
      for jc = 1:numel(cps)
        cp = cps(jc);
        for d = 1:nv
            al = basis_list(c,1+d);     be = basis_list(cp,1+d);
            if kind(d)==3,  ng = 3;     else,   ng = 1;     end
            Td{ic,jc,d} = dir_matrix(kind(d),al,be,PL(c,d),PR(cp,d),JL(c,d),JR(cp,d),...
                                     @(kd,a1,b1,n) getlist(d,kd,a1,b1,n),ELmax(d),ERmax(d),ng);
            [~,cc] = find(Td{ic,jc,d});
            x = cc-1;
            eR = mod(x,ERmax(d)+1);     x = floor(x/(ERmax(d)+1));
            eL = mod(x,ELmax(d)+1);
            hitL{d}(eL+1) = true;       hitR{d}(eR+1) = true;
        end
      end
    end
    ZL = cell(1,numel(ok));     ZR = cell(1,numel(op));
    posL = cell(1,nv);          posR = cell(1,nv);
    for i = 1:numel(ok)
        d = ok(i);  ZL{i} = find(hitL{d})-1;    ZL{i} = ZL{i}(:);
        posL{d} = zeros(ELmax(d)+1,1);  posL{d}(ZL{i}+1) = 1:numel(ZL{i});
    end
    for i = 1:numel(op)
        d = op(i);  ZR{i} = find(hitR{d})-1;    ZR{i} = ZR{i}(:);
        posR{d} = zeros(ERmax(d)+1,1);  posR{d}(ZR{i}+1) = 1:numel(ZR{i});
    end
    nZL = cellfun(@numel,ZL);   nZR = cellfun(@numel,ZR);
    NL = prod([nZL,1]);         NR = prod([nZR,1]);
    nAB = mk(k)*NL*mk(kp)*NR;
    strL = zeros(1,nv);     strR = zeros(1,nv);         % first own variable slowest
    for i = 1:numel(ok),    strL(ok(i)) = prod([nZL(i+1:end),1]);    end
    for i = 1:numel(op),    strR(op(i)) = prod([nZR(i+1:end),1]);    end
    % Pass 2: the triplets of every pair, over the Kronecker product of the
    % per-direction matrices, direction 1 slowest.
    I = cell(ncell,1);  J = cell(ncell,1);  Vv = cell(ncell,1);
    for q = 1:ncell,    I{q} = {};  J{q} = {};  Vv{q} = {};   end
    rsz = zeros(1,nv);  csz = zeros(1,nv);
    for ic = 1:numel(cs)
      c = cs(ic);
      for jc = 1:numel(cps)
        cp = cps(jc);
        T = Td{ic,jc,1};    rsz(1) = size(T,1);    csz(1) = size(T,2);
        for d = 2:nv
            T = kron(T,Td{ic,jc,d});
            rsz(d) = size(Td{ic,jc,d},1);   csz(d) = size(Td{ic,jc,d},2);
        end
        [ri,ci,vv] = find(T);
        if isempty(ri),     continue,   end
        ri = ri(:);     ci = ci(:);     vv = vv(:);
        nt = numel(ri);
        % decode rows to (p,q,j,l) and columns to (gamma,eL,eR) per direction
        p = zeros(nt,nv);   qq = zeros(nt,nv);  j = zeros(nt,nv);   l = zeros(nt,nv);
        gm = zeros(nt,nv);  eL = zeros(nt,nv);  eR = zeros(nt,nv);
        x = ri-1;   y = ci-1;
        for d = nv:-1:1
            rd = mod(x,rsz(d));     x = floor(x/rsz(d));
            l(:,d) = mod(rd,JR(cp,d)+1);    rd = floor(rd/(JR(cp,d)+1));
            j(:,d) = mod(rd,JL(c,d)+1);     rd = floor(rd/(JL(c,d)+1));
            qq(:,d) = mod(rd,PR(cp,d)+1);   p(:,d) = floor(rd/(PR(cp,d)+1));
            cd = mod(y,csz(d));     y = floor(y/csz(d));
            eR(:,d) = mod(cd,ERmax(d)+1);   cd = floor(cd/(ERmax(d)+1));
            eL(:,d) = mod(cd,ELmax(d)+1);   gm(:,d) = floor(cd/(ELmax(d)+1))+1;
        end
        % membership on the (possibly capped) target grids
        t  = tlL{c}(1+[p,j(:,ok)]*stL{c}(:));
        tp = tlR{cp}(1+[qq,l(:,op)]*stR{cp}(:));
        keep = t>0 & tp>0;
        if ~any(keep),  continue,   end
        t = t(keep);    tp = tp(keep);  vv = vv(keep);
        eL = eL(keep,:);    eR = eR(keep,:);    gm = gm(keep,:);
        % coefficient column: (r, left monomial, rp, right monomial)
        lidx = ones(nnz(keep),1);   ridx = ones(nnz(keep),1);
        for i = 1:numel(ok),    d = ok(i);  lidx = lidx+(posL{d}(eL(:,d)+1)-1)*strL(d);    end
        for i = 1:numel(op),    d = op(i);  ridx = ridx+(posR{d}(eR(:,d)+1)-1)*strR(d);    end
        qcell = cell_of_gamma(gm(:,D3));
        % expand over the matrix rows r (space k) and rp (space kp)
        r = 1:mk(k);    rp = reshape(1:mk(kp),1,1,[]);
        pos = (offL(c)+(r-1)*TL(c)+t) + (offR(cp)+(rp-1)*TR(cp)+tp-1)*NLt;
        col = ((r-1)*NL+lidx) + ((rp-1)*NR+ridx-1)*(mk(k)*NL);
        nrep = mk(k)*mk(kp);
        val = repmat(vv,[1,mk(k),mk(kp)]);
        qrep = repmat(qcell,[1,mk(k),mk(kp)]);
        pos = pos(:);   col = col(:);   val = val(:);   qrep = qrep(:);
        if ncell==1
            I{1}{end+1} = pos;  J{1}{end+1} = col;   Vv{1}{end+1} = val;
        else
            [qs,ord] = sort(qrep);
            bnd = [0; find(diff(qs)); numel(qs)];
            for b = 1:numel(bnd)-1
                sel = ord(bnd(b)+1:bnd(b+1));   q = qs(bnd(b)+1);
                I{q}{end+1} = pos(sel);  J{q}{end+1} = col(sel);    Vv{q}{end+1} = val(sel);
            end
        end
        if nrep>1,  continue,   end     % (nrep only documents the expansion)
      end
    end
    Bpos = cell(ncell,1);
    for q = 1:ncell
        if isempty(I{q}),   Bpos{q} = [];    continue,   end
        Bpos{q} = sparse(vertcat(I{q}{:}),vertcat(J{q}{:}),vertcat(Vv{q}{:}),NLt*NRt,nAB);
    end
    blk{k,kp} = struct('Bpos',{Bpos},'ZL',{ZL},'ZR',{ZR},'nAB',nAB,'n3',n3);
  end
end
PM = struct('blk',{blk},'NLt',NLt,'NRt',NRt);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = dir_matrix(kind,al,be,PLc,PRc,JLc,LLc,getlist,ELmax,ERmax,ng)
% One direction of one pair: rows (p,q,j,l), p slowest, columns
% (gamma,eL,eR), gamma slowest, values the table coefficients.
nr = (PLc+1)*(PRc+1)*(JLc+1)*(LLc+1);   nc = ng*(ELmax+1)*(ERmax+1);
[jj,ll] = ndgrid(0:JLc,0:LLc);  jj = reshape(jj,1,[]);   ll = reshape(ll,1,[]);
I = {};     J = {};     V = {};
for p = 0:PLc
    for q = 0:PRc
        L = getlist(kind,al,be,p+q);
        if isempty(L),  continue,   end
        rowidx = ((p*(PRc+1)+q)*(JLc+1)+jj)*(LLc+1)+ll+1;       % 1 x nJL
        eL = L(:,2)+jj;     eR = L(:,3)+ll;                     % nL x nJL
        colidx = ((L(:,1)-1)*(ELmax+1)+eL)*(ERmax+1)+eR+1;
        I{end+1} = repmat(rowidx,size(L,1),1);                              %#ok<AGROW>
        J{end+1} = colidx;                                                  %#ok<AGROW>
        V{end+1} = repmat(L(:,4),1,numel(jj));                              %#ok<AGROW>
    end
end
if isempty(I)
    T = sparse(nr,nc);
    return
end
I = cellfun(@(x) x(:),I,'UniformOutput',false);
J = cellfun(@(x) x(:),J,'UniformOutput',false);
V = cellfun(@(x) x(:),V,'UniformOutput',false);
T = sparse(vertcat(I{:}),vertcat(J{:}),vertcat(V{:}),nr,nc);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ELt,ERt,SE,SF] = target_grids(E,nv,shL,shR)
% The left and right target grids of every component under the shifts,
% with the maps from the source grid.
Nb = numel(E);
ELt = cell(1,Nb);   ERt = cell(1,Nb);   SE = cell(1,Nb);   SF = cell(1,Nb);
for c = 1:Nb
    [ELt{c},SE{c}] = shift_grid(E{c},nv,shL);
    [ERt{c},SF{c}] = shift_grid(E{c},nv,shR);
end
end
