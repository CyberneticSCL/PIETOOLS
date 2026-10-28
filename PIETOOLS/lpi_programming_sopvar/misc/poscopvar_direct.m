function [prog,Pop,Qcell,info] = poscopvar_direct(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,INFO] = POSCOPVAR_DIRECT(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares in the LPI program PROG one positive semidefinite Gram matrix
% Q_c of decision variables per Positivstellensatz term c, and returns the
% self-adjoint 'cdopvar' decision operator on X = L_2^{m_1}[s^1] x ... x
% L_2^{m_M}[s^M] (an empty s^k is R^{m_k})
%
%       Pop = sum_c sum_{b,b'} Z_b'*(g_c*Q_c{b,b'})*Z_{b'} >= 0,
%
% with Z_b the basis operators of all spaces at the degrees of term c and
% g_c(theta) >= 0 the weight of term c. Basis operator b of space k
% carries a multi-index alpha over the variables of s^k and maps
% L_2^{m_k}[s^k] to L_2[theta], theta over every registry variable:
%
%       (Z_b x)(theta) = int I_alpha(theta-s) ...
%                            (I_{m_k} kron Z^alpha(theta,s)) x(s) ds,
%
% I_alpha(r) = prod_d I_{alpha_d}(r_d) over the variables of s^k, alpha_d =
% 1 the multiplier (delta), 2 the lower integral (s_d <= theta_d), 3 the
% upper integral (s_d >= theta_d), 4 the full integral; Z^alpha(theta,s)
% is the vector of monomials of the degrees DEG, and on R^{m_k},
% Z_b x = (I_{m_k} kron Z^alpha(theta)) x. The kernel coefficients of Pop
% are mapped from the entries of each Q_c directly: with Gram row
% u = (b, matrix row r, monomial theta^p s^j), column v = (b', t,
% theta^q s'^l), alpha and beta the multi-indices of b and b', and
% g_c = prod_d g_d, the kernel of Pop is the sum over (u,v) of
%
%   Q_uv s^j s'^l e_r e_t' prod_d int I_{alpha_d}(theta_d-s_d) ...
%        I_{beta_d}(theta_d-s'_d) g_d(theta_d) theta_d^(p_d+q_d) dtheta_d.
%
% POSCOPVAR_DIRECT('clear') empties the cache of position maps and returns
% [], [], {} and struct('cleared',true).
%
% INPUT
% - prog:     'struct', the LPI program (LPIPROGRAM_SOP), or the char
%             'clear';
% - dims:     double, M x 1 positive integers, dims(k) = m_k the number of
%             components of space k; a scalar is used for every space;
% - spaces:   1 x M cell of 'cellstr', spaces{k} the names of the
%             variables of s^k ({} for R^{m_k}); a plain 'cellstr' is one
%             space;
% - dom:      double, nv x 2, dom(d,:) the interval of variable d of the
%             registry (the sorted union of all spaces), or one 1 x 2 row
%             for every variable; or a 'struct' with fields vars
%             ('cellstr') and dom (double, one row per name);
% - deg:      the degrees of Z^alpha(theta,s): a scalar (int and mult), or
%             a 'struct' with the optional fields int (scalar or 1 x nv,
%             the degree in each theta_d, default 1; alias 'weight'), mult
%             (scalar or 1 x nv, the degree in each s_d of s^k, 0 where
%             alpha_d = 1, default int; alias 'lift'), joint (scalar, a cap
%             on the total degree) and subset (2^(2nv) caps, one per
%             subset of [theta, s^k], for a space with every registry
%             variable only); or a 1 x M cell of these, one per space, each
%             of which may be a cell with one per basis operator of the
%             space;
% - options:  (optional) 'struct' with the optional fields
%               type          char, 'pos' (any other value is an error);
%               sep           logical, scalar or 1 x nv: where true, the
%                             lower and upper basis operators of that
%                             direction are replaced by the full integral
%                             (alpha_d = 4);
%               include       1 x M cell, entry k selecting basis
%                             operators of space k: multi-index rows over
%                             the variables of s^k, a logical mask, or
%                             linear indices (default all);
%               psatz         double, the term codes, one or a row, one
%                             Gram each in the given order: 0 (default,
%                             g = 1), 1 (g = prod_d (theta_d-a_d)*
%                             (b_d-theta_d)), 2d+1 / 2d+2
%                             (g = (theta_d-a_d)/L_d or (b_d-theta_d)/L_d,
%                             L_d = b_d-a_d), 2nv+2+d
%                             (g = (theta_d-a_d)*(b_d-theta_d)), d = 1..nv;
%               psatz_offset  double, nonnegative integers, one per code or
%                             one for all (default 0): the int of the term
%                             lowered by it in every direction (in
%                             direction d alone for code 2nv+2+d), floored
%                             at 0, mult kept, a joint cap lowered by it
%                             per direction the weight acts in; an error
%                             with a subset cap;
%               path          char, 'fast' (default) or 'general': the
%                             weight moved onto the Gram index as a shift
%                             when the map of the shifted index set is
%                             cached, else integrated with the monomials;
%                             the same Pop up to rounding;
%               cache         logical (default true): keep the position
%                             maps for the session, keyed on everything
%                             they depend on;
%               estimate      logical (default false): return the sizes
%                             only, declaring nothing.
%
% OUTPUT
% - prog:     'struct', the program with each Q_c declared positive
%             semidefinite (SOSQUADVAR); unchanged with options.estimate;
% - Pop:      'cdopvar' object, M x M blocks on X, the operator above; []
%             with options.estimate;
% - Qcell:    N_b x N_b cell, Qcell{b,b'} the 'cellstr' names of the
%             decision variables of block (b,b') of Q_c, N_b the number of
%             basis operators; with several terms, a 1 x nt cell of such
%             cells, one per term; {} with options.estimate;
% - info:     'struct' with fields
%               basis_list    double, N_b x (1+nv), row b: the space of
%                             basis operator b, then its multi-index entry
%                             alpha_d in each registry direction, 0 where
%                             its space lacks variable d;
%               E             1 x N_b cell, the exponents of
%                             Z^alpha(theta,s) of each basis operator,
%                             first term;
%               N             double row, the Gram dimension of each term;
%               ndec          double, the number of decision variables of
%                             Pop;
%               path          char, the path used;
%               route         char per term: 'fast-shift' (the shift on a
%                             cached map), 'fast-own' (the map of the term
%                             itself) or 'general'; a cell with several
%                             terms;
%               code, offset  double rows, the terms;
%               cache_hit     logical row, per term;
%               Qdpvar        'dpvar' object, the N x N Gram matrix, for
%                             SOSGETSOL; a cell with several terms;
%               terms         struct array, per term: code, offset, E, N,
%                             off, Nc, names, ndec, route, cache_hit,
%                             t_declare, t_map, t_product;
%               t_declare, t_map, t_product, t_blocks, t_total
%                             double, seconds;
%             with options.estimate, the fields N, ndec (sum(N.*(N+1)/2))
%             and estimate (true).
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
% MMP, 10/08/2026: Single-direction box quadratics as term codes 2nv+2+d,
%                offset lowering 'int' in that direction only, so the
%                Markov-Lukacs tensor set (plain, g_d at w-1 in s_d, the
%                product at w-1) can be compared with the 2N faces at w.
% MMP, 10/08/2026: Several Positivstellensatz terms in one call ('psatz' a
%                row, 'psatz_offset' per term): one Gram per term, the maps
%                shared where the index sets coincide, every term's map
%                embedded in the union of the bases and summed on one
%                decision list, so a sum of terms costs one block build and
%                no container addition. For 'lpi_ineq_sop'. One term is
%                unchanged (the same program, Qcell and info as before).
% MMP, 10/17/2026: Spaces without a spatial variable (every space R^q, nv =
%                0: the W on R^{n_w} or R^{n_z} of the H2 executives through
%                poslpivar_sop -> poscopvar) errored in pos_map, found by
%                the Posopvar session's chain test: the stride vector of a
%                one-row exponent grid had a stray entry, and the Kronecker
%                product over directions indexed a third dimension of size
%                zero. Both now read the empty direction set as one Gram
%                position mapped to one coefficient with weight 1; the map
%                is then the Gram itself, block by block, as copquadvar
%                builds it. Nothing changes when nv >= 1.
% MMP, 10/21/2026: options.estimate: the Gram dimension of every term and
%                the decision variables of the call, from the grids the
%                build uses, returned before any declaration. LPI_INEQ_SOP
%                prints them before it builds (maintainer: lpi_ineq_sop
%                always estimates the decision variables needed).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body and headers
%                of the local functions; code unchanged. The deleted
%                header text was the comparison with copquadvar, the
%                references to poscopvar and copquadvar for the inputs,
%                the reason for the cache rule of the fast path, a test
%                name and the cost note; the 1-D table of the map is in
%                the header of dir_list. The help text of the 10/21/2026
%                entry is in the entries of options.estimate and of the
%                outputs, without markers.

% The cache of position maps; 'clear' empties it
persistent CACHE
if nargin>=1 && (ischar(prog) || isstring(prog))
    if strcmpi(char(prog),'clear')
        CACHE = [];     prog = [];  Pop = [];  Qcell = {};  info = struct('cleared',true);
        return
    end
    error("The first input is the program, or 'clear' to empty the cache.")
end
% Check the inputs; OPTIONS a struct, of type 'pos' only
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
% The terms: codes and nonnegative integer offsets, one per code
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
% Path, cache and estimate options; the timer
path = 'fast';
if isfield(options,'path') && ~isempty(options.path),    path = lower(char(options.path));   end
if ~any(strcmp(path,{'fast','general'}))
    error("'path' should be 'fast' or 'general'.")
end
usecache = true;
if isfield(options,'cache') && ~isempty(options.cache),  usecache = logical(options.cache);  end
est = isfield(options,'estimate') && ~isempty(options.estimate) && logical(options.estimate); % MMP, 10/21/2026
Nest = zeros(1,nt);             % the Gram dimensions, when est             % MMP, 10/21/2026
t_all = tic;

% % % Spaces, registry, components: 'copquadvar's enumeration, through
% 'lift_copvar' (multiindex_grid, direction 1 fastest, then include and sep).
[meta,sp] = parse_copvar_spaces(dims,spaces,dom);
M = numel(sp);  vars = meta.vars;   nv = numel(vars);   dom = meta.dom;
mask = meta.space_out;  mk = meta.dim_out(:);
own = cell(1,M);
for k = 1:M,    own{k} = reshape(find(mask(k,:)),1,[]);     end
% Term codes valid for nv variables
% if any(~ismember(codes,[0,1,3:2*nv+2]))                                   % MMP, 10/08/2026 (was)
if any(~ismember(codes,[0,1,3:3*nv+2]))                                     % MMP, 10/08/2026
    error("'psatz' entries should be 0, 1, face codes 2d+1 / 2d+2 or quadratic codes 2nv+2+d with d in 1..nv.")
end
% The basis operators of every space, their space and index within it
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
%       spec = reduce_spec(spec,offs(c),ng_of(codes(c),nv),nv);             % MMP, 10/08/2026 (was)
        offv = offs(c);                                                     % MMP, 10/08/2026
        if codes(c)>2*nv+2      % single-direction quadratic: lower 'int' there only
            offv = zeros(1,nv);     offv(codes(c)-2*nv-2) = offs(c);        % MMP, 10/08/2026
        end
        spec = reduce_spec(spec,offv,ng_of(codes(c),nv),nv);                % MMP, 10/08/2026
        % Caps per variable, 'mult' zero where the multi-index is 1; the
        % exponent grid
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
    if est,     Nest(c) = N;    continue,   end     % sizes only            % MMP, 10/21/2026
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
    % The route taken
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
    % Record the term
    terms(c) = struct('code',codes(c),'offset',offs(c),'E',{E},'N',N,'off',off,'Nc',Nc,...
                      'names',{names},'ndec',ndec_c,'route',route,'cache_hit',hit,...
                      'Qdpvar',Pq{1,1},'PM',PM,'SE',{SE},'SF',{SF},'useshift',useshift,...
                      't_declare',t_declare,'t_map',t_map,'t_product',0);
end
if est                          % one Gram per term, distinct names         % MMP, 10/21/2026
    Pop = [];   Qcell = {};                                                 % MMP, 10/21/2026
    info = struct('N',Nest,'ndec',sum(Nest.*(Nest+1)/2),'estimate',true);   % MMP, 10/21/2026
    return                                                                  % MMP, 10/21/2026
end                                                                         % MMP, 10/21/2026

% % % The decision variables of all terms in one list, sorted as copquadvar
% sorts its names, and each term's placement: dvar k sits at (u,v) and
% (v,u) of its Gram, rows in the common list; the shifts move the source
% positions to the target positions.
allnames = vertcat(terms.names);
[Zd,~,rowall] = unique(allnames);
ndec = numel(Zd);
D = cell(1,nt);     pos = 0;
for c = 1:nt
    % Rows of the names of the term in the common list; the name of each
    % Gram position (u,v)
    t0 = tic;
    rm = rowall(pos+(1:terms(c).ndec));     pos = pos+terms(c).ndec;
    N = terms(c).N;
    [U,V] = ndgrid(1:N,1:N);
    kk = lind_sym(U(:),V(:),N);
    D{c} = sparse(rm(kk),U(:)+(V(:)-1)*N,1,ndec,N^2);
    % Shift route: from the source positions to the target positions
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
        % Block (k,kp): the union of the bases the terms reach
        ok = own{k};    op = own{kp};
        ZL = cell(1,numel(ok));     ZR = cell(1,numel(op));
        for c = 1:nt
            b = terms(c).PM.blk{k,kp};
            for i = 1:numel(ok),    ZL{i} = union(ZL{i},b.ZL{i});   ZL{i} = ZL{i}(:);  end
            for i = 1:numel(op),    ZR{i} = union(ZR{i},b.ZR{i});   ZR{i} = ZR{i}(:);  end
        end
        % Its coefficient count and parameter cells
        NL = prod([cellfun(@numel,ZL),1]);  NR = prod([cellfun(@numel,ZR),1]);
        nAB = mk(k)*NL*mk(kp)*NR;
        n3 = terms(1).PM.blk{k,kp}.n3;      ncell = 3^n3;
        % Per cell, the sum over the terms of names times positions, on
        % the union bases
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
        % Parameters: a zero constant part, B the decision coefficients
        params = struct('A',{cell(ncell,1)},'B',{cell(ncell,1)});
        for q = 1:ncell
            params.A{q} = sparse(nAB,1);
            if isempty(Bsum{q}),    params.B{q} = sparse(ndec,nAB);
            else,                   params.B{q} = Bsum{q};
            end
        end
        params.A = reshape(params.A,[3*ones(1,n3),1,1]);
        params.B = reshape(params.B,[3*ones(1,n3),1,1]);
        % The 'sdopvar' block on the variables and domains of the spaces
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
% Route and Gram per term; the info output
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
% N = NG_OF(CODE,NV) returns the number of directions in which the weight
% of the term CODE has positive degree: 0 for the plain term, NV for the
% box product, 1 for a face or a single-direction box quadratic.
% INPUT  code: double, the term code; nv: double, the number of registry
%        variables.
% OUTPUT n:    double.
if code==0,         n = 0;
elseif code==1,     n = nv;
else,               n = 1;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spec = reduce_spec(spec,off,ng,nv)
% SPEC = REDUCE_SPEC(SPEC,OFF,NG,NV) returns the degree specification of a
% Psatz term: 'int' lowered by OFF in every direction (OFF a row: per
% direction), floored at 0, 'mult' kept at what it was (its default is the
% ORIGINAL 'int'), a 'joint' cap lowered by OFF per direction the weight
% acts in (a row: by its sum). A 'subset' cap cannot move with 'int' and
% is refused.
% INPUT  spec: double scalar, or 'struct' with the optional fields int,
%              mult, joint, subset;
%        off:  double, the offset: a scalar, or 1 x nv, one per direction;
%        ng:   double, the number of directions the weight acts in;
%        nv:   double, the number of registry variables.
% OUTPUT spec: SPEC when OFF is zero, else a 'struct' with fields int
%              (1 x max(nv,1) double), mult and, if given, joint.
% No offset: SPEC unchanged
% if off==0,  return,     end                                               % MMP, 10/08/2026 (was)
if all(off==0),     return,     end                                         % MMP, 10/08/2026
% A scalar specification as a struct; a subset cap refused
if isnumeric(spec) && isscalar(spec)
    spec = struct('int',spec,'mult',spec);
elseif ~isstruct(spec)
    error("Monomial degrees should be specified as a scalar or 'struct' object.")
end
if isfield(spec,'subset') && ~isempty(spec.subset)
    error("A 'psatz_offset' is not defined for a 'subset' degree specification.")
end
% 'int' lowered and floored at 0; 'mult' defaults to the original 'int'
w = 1;
if isfield(spec,'int') && ~isempty(spec.int),   w = spec.int;   end
if ~isfield(spec,'mult') || isempty(spec.mult), spec.mult = w;  end
w = reshape(w,1,[]);    if isscalar(w),   w = repmat(w,1,max(nv,1));    end
spec.int = max(w-off,0);
% 'joint' lowered by the total offset over the directions of the weight
if isfield(spec,'joint') && ~isempty(spec.joint)
%   spec.joint = max(spec.joint-off*ng,0);                                  % MMP, 10/08/2026 (was)
    tot = off*ng;   if ~isscalar(off),  tot = sum(off);     end             % MMP, 10/08/2026
    spec.joint = max(spec.joint-tot,0);                                     % MMP, 10/08/2026
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Ec = embed_cols(ZLc,ZRc,ZL,ZR,m,n)
% EC = EMBED_COLS(ZLC,ZRC,ZL,ZR,M,N) returns the sparse map from the
% coefficient columns of a term, on its own bases ZLc, ZRc, to the columns
% on the union bases ZL, ZR of the block; [] when the bases coincide.
% Columns are (row r, left monomial, column t, right monomial), left and
% right monomials first variable slowest.
% INPUT  ZLc, ZRc: cells of exponent columns, one per output / input
%                  variable of the block: the bases of the term;
%        ZL, ZR:   cells of the same lengths, the union bases;
%        m, n:     double, the component counts of the output and input
%                  spaces of the block.
% OUTPUT Ec:       sparse 0/1 matrix, m*NLc*n*NRc x m*NL*n*NR, or [].
% Same bases: no map
if isequal(ZLc,ZL) && isequal(ZRc,ZR),  Ec = [];  return,     end
% The position of each own monomial in the union bases
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
% The union column of every own column
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
% K = LIND_SYM(U,V,N) returns the decision variable of Gram position
% (u,v): the numbering of SOSQUADVAR of the lower triangle by columns,
% mirrored.
% INPUT  u, v: double arrays of one size, row and column indices in 1..N;
%        N:    double, the Gram dimension.
% OUTPUT k:    double array of that size, indices in 1..N*(N+1)/2.
U = max(u,v);   V = min(u,v);
k = U + (V-1)*N - V.*(V-1)/2;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = blkdiag_comp(Sc,m)
% S = BLKDIAG_COMP(SC,M) returns the shift of the whole Gram index from the
% per-component shifts: rows (component, matrix row, monomial), so
% component c contributes kron(I_m, Sc{c}).
% INPUT  Sc: 1 x Nb cell of sparse matrices, the shift map of each basis
%            operator;
%        m:  Nb x 1 double, the component count of the space of each.
% OUTPUT S:  sparse, blkdiag of kron(speye(m(c)),Sc{c}) over c.
B = cell(1,numel(Sc));
for c = 1:numel(Sc),    B{c} = kron(speye(m(c)),Sc{c});   end
S = blkdiag(B{:});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg = alias_deg(deg)
% DEG = ALIAS_DEG(DEG) returns the degree specification DEG with the fields
% 'lift' and 'weight' (the names of POSCOPVAR_LIFT) read as 'mult' and
% 'int': each fills the field it stands for when that is absent or empty,
% and is removed.
% INPUT  deg: double, 'struct', or a cell of them, nested to any depth.
% OUTPUT deg: the same, without the fields 'lift' and 'weight'.
% A cell entry by entry; a struct field by field
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
% [GCO,SHL,SHR] = TERM_WEIGHT(CODE,NV,DOM,PATH) returns the weight of the
% term CODE per direction: gco{d}, the coefficients of g_d (constant
% first) for the general path (1 on the fast path); shL{d}, shR{d}, the
% rows [shift, coefficient] of the left and right degree-1 factors on the
% fast path (the identity [0,1] on the general path).
% INPUT  code: double, the term code (0, 1, 2d+1, 2d+2 or 2nv+2+d);
%        nv:   double, the number of registry variables;
%        dom:  nv x 2 double, the registry domains;
%        path: char, 'fast' or 'general'.
% OUTPUT gco:      1 x nv cell of double rows;
%        shL, shR: 1 x nv cells of K x 2 double.
% The identity weight; then the factors of CODE
gco = repmat({1},1,nv);     shL = repmat({[0,1]},1,nv);   shR = shL;
if code==1
    % Box product in every direction: on the fast path (theta_d-a_d) on the
    % left and (b_d-theta_d) on the right, else their product
    for d = 1:nv
        a = dom(d,1);   b = dom(d,2);
        if strcmp(path,'fast')
            shL{d} = [0,-a; 1,1];       shR{d} = [0,b; 1,-1];
        else
            gco{d} = [-a*b, a+b, -1];
        end
    end
elseif code>2*nv+2      % the box quadratic in direction d alone            % MMP, 10/08/2026
    d = code-2*nv-2;    a = dom(d,1);   b = dom(d,2);                       % MMP, 10/08/2026
    if strcmp(path,'fast'),     shL{d} = [0,-a; 1,1];   shR{d} = [0,b; 1,-1]; % MMP, 10/08/2026
    else,                       gco{d} = [-a*b, a+b, -1];                   % MMP, 10/08/2026
    end
elseif code>=3
    % Face of direction d: (theta_d-a_d)/L or (b_d-theta_d)/L, on the fast
    % path a shift on the left only
    d = floor((code-1)/2);  a = dom(d,1);   b = dom(d,2);   L = b-a;
    if mod(code,2)==1,  fl = [-a,1]/L;  else,   fl = [b,-1]/L;    end
    if strcmp(path,'fast'),     shL{d} = [0,fl(1); 1,fl(2)];
    else,                       gco{d} = fl;
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Et,S] = shift_grid(E,nv,sh)
% [ET,S] = SHIFT_GRID(E,NV,SH) returns the target grid of a component under
% the per-direction shifts sh{d} (rows [e, coef]) and the sparse map S,
% S(t_source,t_target) = the product of the coefficients of the shift
% combination taking one to the other.
% INPUT  E:   T x (nv+n_k) double, the exponents of a basis operator, the
%             theta exponents first;
%        nv:  double, the number of registry variables;
%        sh:  1 x nv cell of K_d x 2 double, the shifts.
% OUTPUT Et:  double, the distinct shifted rows, sorted;
%        S:   sparse T x size(Et,1).
% Every combination of one shift per direction
T = size(E,1);
combos = multiindex_grid(cellfun(@(s) 1:size(s,1),sh,'UniformOutput',false));
nc = size(combos,1);
% The shifted rows, their source rows and coefficients
Eall = zeros(nc*T,size(E,2));   rows = zeros(nc*T,1);    vals = zeros(nc*T,1);
for i = 1:nc
    ev = zeros(1,nv);   cf = 1;
    for d = 1:nv
        ev(d) = sh{d}(combos(i,d),1);   cf = cf*sh{d}(combos(i,d),2);
    end
    Es = E;     Es(:,1:nv) = Es(:,1:nv)+ev;
    Eall((i-1)*T+(1:T),:) = Es;     rows((i-1)*T+(1:T)) = (1:T)';   vals((i-1)*T+(1:T)) = cf;
end
% The distinct targets and the map
[Et,~,ic] = unique(Eall,'rows');
S = sparse(rows,ic,vals,T,size(Et,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key = cache_key(basis_list,ELt,ERt,mk,dom,gco,own)
% KEY = CACHE_KEY(BASIS_LIST,ELT,ERT,MK,DOM,GCO,OWN) returns the key of a
% position map in the cache: the JSON text of everything the map depends
% on.
% INPUT  basis_list: Nb x (1+nv) double, the basis operators;
%        ELt, ERt:   1 x Nb cells, the left and right target grids;
%        mk:         M x 1 double, the component counts of the spaces;
%        dom:        nv x 2 double, the registry domains;
%        gco:        1 x nv cell, the coefficients of the weight;
%        own:        1 x M cell, the registry indices of each space.
% OUTPUT key:        char.
s = jsonencode(struct('bl',basis_list,'EL',{ELt},'ER',{ERt},'m',mk,'dom',dom,...
                      'g',{gco},'own',{own}));
key = s;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = dir_list(kind,al,be,n,g,a,b)
% L = DIR_LIST(KIND,AL,BE,N,G,A,B) returns, for one direction and the
% multi-index entries AL (left) and BE (right) of a pair of basis
% operators, int I_AL(theta-s) I_BE(theta-s') G(theta) dtheta with
% G = g(theta) theta^N, as rows [gamma, left exponent, right exponent,
% coefficient]: the kernel, coefficient times s^left s'^right, in
% parameter cell gamma (1 multiplier, 2 lower s' <= s, 3 upper s <= s';
% 1 when the direction is not shared). With calG the antiderivative of G,
%
%   (1,1) mult: G(s)        (1,2) low: G(s)            (1,3) up: G(s)
%   (2,1) up:   G(s')       (3,1) low: G(s')
%   (2,2) low:  calG(b)-calG(s)      up: calG(b)-calG(s')
%   (3,3) low:  calG(s')-calG(a)     up: calG(s)-calG(a)
%   (2,3) up:   calG(s')-calG(s)     (3,2) low: calG(s)-calG(s')
%
% an index 4 (full) is the sum of the rows (or columns) 2 and 3, and a
% direction only one space has, or neither has, gives the one-sided or
% the definite integral of G.
% INPUT  kind:   double, 3 the direction is in both spaces, 2 in the
%                output space only, 1 in the input space only, 0 in
%                neither;
%        al, be: double, 1 to 4, the multi-index entries of the left and
%                right basis operators in the direction;
%        n:      double, the exponent of theta;
%        g:      double row, the coefficients of g, constant first;
%        a, b:   double, the domain of the direction.
% OUTPUT L:      K x 4 double, the rows with a nonzero coefficient.
% G(s), G(s'), calG(s), calG(s'), calG(b) and calG(a) as rows
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
% The entry of the table for KIND, AL and BE
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
% Only the nonzero coefficients
L = L(L(:,4)~=0,:);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PM = pos_map(ELt,ERt,basis_list,mk,own,dom,gco,nv)
% PM = POS_MAP(ELT,ERT,BASIS_LIST,MK,OWN,DOM,GCO,NV) returns the map from
% the target Gram positions to the coefficient columns of every block and
% parameter cell: PM.blk{k,kp} with fields Bpos (cell over the 3^n3 cells,
% each sparse NLt*NRt x nAB), ZL, ZR, nAB, n3.
% INPUT  ELt, ERt:   1 x Nb cells, the left and right target grids of the
%                    basis operators (exponent rows over [theta, s^k]);
%        basis_list: Nb x (1+nv) double, the basis operators;
%        mk:         M x 1 double, the component counts of the spaces;
%        own:        1 x M cell, the registry indices of each space;
%        dom:        nv x 2 double, the registry domains;
%        gco:        1 x nv cell, the coefficients of g_d, constant first;
%        nv:         double, the number of registry variables.
% OUTPUT PM:         'struct' with fields blk (M x M cell of 'struct':
%                    Bpos, a 3^n3 x 1 cell of sparse matrices or []; ZL,
%                    ZR, the monomial bases reached; nAB, the number of
%                    coefficients; n3, the number of shared variables),
%                    NLt and NRt (the numbers of left and right Gram
%                    positions).
% Sizes and offsets of the Gram rows and columns of each basis operator
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
    % monomial index on the target grid from the exponent row, 0 if absent.
    % The stride vector has one entry per exponent column: with no spatial
    % variable (every space R^q) it is EMPTY, and the one-row grid maps to
    % index 1; cumprod([1,rng(1:end-1)]) gave a stray entry there.         % MMP, 10/17/2026
%   rng = [PL(c,:)+1, JL(c,ok)+1];      stL{c} = cumprod([1,rng(1:end-1)]); % MMP, 10/17/2026 (was)
    rng = [PL(c,:)+1, JL(c,ok)+1];      st = cumprod([1,rng]);  stL{c} = st(1:numel(rng)); % MMP, 10/17/2026
    tlL{c} = zeros(prod(rng),1);        tlL{c}(1+ELt{c}*stL{c}(:)) = (1:TL(c))';
%   rng = [PR(c,:)+1, JR(c,ok)+1];      stR{c} = cumprod([1,rng(1:end-1)]); % MMP, 10/17/2026 (was)
    rng = [PR(c,:)+1, JR(c,ok)+1];      st = cumprod([1,rng]);  stR{c} = st(1:numel(rng)); % MMP, 10/17/2026
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
        % L = GETLIST(D,KIND,AL,BE,N) returns DIR_LIST for direction D,
        % memoized in LST on (KIND,AL,BE,N).
        % INPUT  d: double, the direction; kind, al, be, n: double, as
        %        DIR_LIST takes them.
        % OUTPUT L: K x 4 double (zeros(0,4) when empty).
        L = LST{d}{kind+1,al+1,be+1,n+1};
        if isempty(L)
            L = dir_list(kind,al,be,n,gco{d},dom(d,1),dom(d,2));
            if isempty(L),  L = zeros(0,4);     end
            LST{d}{kind+1,al+1,be+1,n+1} = L;
        end
    end

% Every block (k,kp): the kind of each direction (3 shared, 2 output
% space only, 1 input space only, 0 neither) and its basis operators
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
    % The bases: the exponents hit, and the position of each in them
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
    % The number of coefficients; strides of the left and right monomials
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
        % With no spatial variable the Kronecker product over directions is
        % empty: one Gram position (the constant) maps to one coefficient,
        % weight 1 (Td has no third dimension to index then).              % MMP, 10/17/2026
%       T = Td{ic,jc,1};    rsz(1) = size(T,1);    csz(1) = size(T,2);      % MMP, 10/17/2026 (was)
        if nv==0,   T = sparse(1);                                          % MMP, 10/17/2026
        else,       T = Td{ic,jc,1};    rsz(1) = size(T,1);    csz(1) = size(T,2); % MMP, 10/17/2026
        end                                                                 % MMP, 10/17/2026
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
        % The triplets, grouped by parameter cell
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
    % One sparse matrix per parameter cell
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
% T = DIR_MATRIX(KIND,AL,BE,PLC,PRC,JLC,LLC,GETLIST,ELMAX,ERMAX,NG) returns
% one direction of one pair of basis operators: rows (p,q,j,l), p slowest,
% the theta exponents p, q and the lift exponents j, l of the left and
% right Gram monomials; columns (gamma,eL,eR), gamma slowest, the parameter
% cell and the exponents of s and s'; values the table coefficients.
% INPUT  kind, al, be: double, as DIR_LIST takes them;
%        PLc, PRc:     double, the largest left and right theta exponents;
%        JLc, LLc:     double, the largest left and right lift exponents;
%        getlist:      function handle, L = getlist(kind,al,be,n), the
%                      rows of DIR_LIST;
%        ELmax, ERmax: double, the largest exponents of s and s';
%        ng:           double, the number of parameter cells (3 shared, 1
%                      otherwise).
% OUTPUT T:            sparse (PLc+1)*(PRc+1)*(JLc+1)*(LLc+1) x
%                      ng*(ELmax+1)*(ERmax+1).
% Sizes; the lift exponents (j,l) of every pair
nr = (PLc+1)*(PRc+1)*(JLc+1)*(LLc+1);   nc = ng*(ELmax+1)*(ERmax+1);
[jj,ll] = ndgrid(0:JLc,0:LLc);  jj = reshape(jj,1,[]);   ll = reshape(ll,1,[]);
% Triplets for every theta pair (p,q): the table at exponent p+q, shifted
% by the lift exponents
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
% No triplet: the zero matrix; else one sparse matrix
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
% [ELT,ERT,SE,SF] = TARGET_GRIDS(E,NV,SHL,SHR) returns the left and right
% target grids of every component under the shifts, with the maps from
% the source grid (SHIFT_GRID).
% INPUT  E:        1 x Nb cell, the exponent grid of each basis operator;
%        nv:       double, the number of registry variables;
%        shL, shR: 1 x nv cells of K x 2 double, the left and right shifts.
% OUTPUT ELt, ERt: 1 x Nb cells, the target grids;
%        SE, SF:   1 x Nb cells of sparse maps, source rows to target rows.
Nb = numel(E);
ELt = cell(1,Nb);   ERt = cell(1,Nb);   SE = cell(1,Nb);   SF = cell(1,Nb);
for c = 1:Nb
    [ELt{c},SE{c}] = shift_grid(E{c},nv,shL);
    [ERt{c},SF{c}] = shift_grid(E{c},nv,shR);
end
end
