function [prog,Pop,info] = lpi_ineq_sop(prog,P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,INFO] = LPI_INEQ_SOP(PROG,P,OPTS) adds to the LPI program PROG
% the constraint P >= 0 for a self-adjoint PI operator P. For a 'dopvar'
% or 'opvar' P on X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M] (an empty s^k
% is R^{m_k}), it declares the positive operator
%
%       Pop = sum_c Z_c'*(g_c*Q_c)*Z_c,     Q_c >= 0,
%
% one term c per weight g_c(theta) >= 0 on the domain, Z_c a column of
% multiplier and integral operators from X to L_2[theta] with monomial
% degrees read from the support of P or of OPTS.like by the lift/weight
% rule (below), Q_c a Gram matrix of decision variables, and imposes
% P - Pop == 0, one block of each adjoint pair. Before declaring Pop, it
% prints the number of decision variables Pop adds and the program total.
% For a 'dpvar' or 'polynomial' P it imposes P >= 0 as an SOS inequality;
% P of any other class is an error.
%
% The lift/weight rule: in each direction d, with Dmin_d and Dmax_d the
% largest min(i,j) and max(i,j) over the exponents i, j of s_d, s'_d in
% the lower kernels of the support, Mdeg_d the largest multiplier degree
% and J_d the largest degree in a block between R^q and L2, the lift
% degree is D_d = Dmin_d + dD and the weight degree is
%
%   w_d = max(ceil((Dmax_d-D_d-1)/2), ceil(Mdeg_d/2), J_d-D_d-1, 0) + dw,
%
% at least wmin when the support is read from P itself.
%
% INPUT
% - prog:  'struct', the LPI program; read are decvartable and, if
%          present, settings.lpi.ineq, sop_pos (cell of the positive
%          'dopvar' operators declared in PROG) and sopeq (cell of the
%          records of equality rows);
% - P:     'dopvar' or 'opvar' object, square (the same output and input
%          spaces) and self-adjoint (not checked), the operator to make
%          positive; for an 'opvar' P the constraint certifies P >= 0. A
%          'dpvar' or 'polynomial' P is imposed as an SOS inequality; any
%          other class is an error;
% - opts:  (optional) 'struct' with the optional fields below, or an LPI
%          settings structure, of which only lpi.ineq is read. dD, dw and
%          wmin not given are those of the settings of the program,
%          prog.settings.lpi.ineq (the 'light' tier for a program without
%          them), where a vector holds the entries for 1, 2, and 3 or more
%          variables:
%            dD, dw        double, added to the lift degree and to the
%                          weight degree read from the support (default
%                          1, 1);
%            wmin          double, the least weight degree when the
%                          support is read from P itself, not from
%                          opts.like (default 2 in one and two variables,
%                          0 in three or more; 0 when opts.dw is given and
%                          opts.wmin is not);
%            psatz         the terms of Pop: 'auto' (default; in one
%                          variable the plain term and the box product one
%                          weight degree lower, in nv ~= 1 variables the
%                          plain term and the 2nv face terms (s_d-a_d),
%                          (b_d-s_d), one degree lower when every weight
%                          degree is at least 2), 'none' (the plain term),
%                          'product' (as 'auto' in one variable), 'faces'
%                          (the plain term and the 2nv face terms),
%                          'corners' (the face terms and every face
%                          product of two directions), 'faces1' (the face
%                          terms one weight degree lower for every w; in
%                          one variable as 'product'), or a double row of
%                          term codes: 0 plain term, 1 box product, 2d+1 /
%                          2d+2 face (s_d-a_d) / (b_d-s_d), 2nv+2+d box
%                          quadratic of direction d;
%            psatz_offset  double, with a row of codes: the weight degree
%                          by which each term is lowered, one per code or
%                          one for all (default 0);
%            psatz_rows, psatz_rows_offset, psatz_rows_perdir
%                          with a row of codes: R x nv rows of
%                          per-direction factor codes (0 none, 1 face
%                          (s_d-a_d), 2 face (b_d-s_d), 3 box quadratic),
%                          one term per row with the product of its
%                          factors as weight (face products), their
%                          offsets, and logical perdir (a row lowers the
%                          degree only in the directions of its nonzero
%                          factors);
%                          Default for psatz and psatz_offset: the fields
%                          psatz, psatz_offset of the settings tier
%                          (st.lpi.ineq) when present, else 'auto'; for
%                          the three row fields, when the call gives no
%                          psatz: those of the tier, when present;
%                          a tier psatz may be a cell {1-D, 2-D, N-D}, its
%                          entry for the number of registry variables of P;
%            wcross        double, the weight degree of an L2 space in every
%                          registry direction it does not own, in place of
%                          w_d (0: no slack monomial there);
%            dw_space      1 x M double, added to the weight degree of
%                          space k in every direction;
%                          default for both: those of the tier, if any;
%            Dmax, cap_inc logical, the lift degree max(Dmax_d, D_d); and
%                          'auto' or 1 x nv double, the cap theta_d + s_d
%                          <= cap_inc_d on the lower and upper basis
%                          operators ('auto': Dmax_d - 1); default: those
%                          of the tier, if any (false, none);
%            like          'dopvar' or 'opvar' object, or a cell of them,
%                          with variables among those of P: their pooled
%                          support is read instead of that of P. Default
%                          for a 'dopvar' P without deg: the operators of
%                          prog.sop_pos with variables among those of P
%                          whose decision variables P carries (tested on
%                          up to five of them); if there is none, P;
%            deg           a degree specification used instead of the one
%                          read from the support: a scalar, a 'struct' with
%                          fields int, mult, joint, subset, or a 1 x M cell
%                          of them, one per space; the terms are then psatz
%                          and psatz_offset for a numeric psatz, else those
%                          that psatz names;
%            prune         logical (default true): on each space with
%                          spatial variables, keep only the basis operators
%                          whose diagonal term reaches a nonzero parameter
%                          cell of the diagonal block of P;
%            margin        double >= 0 (default 0): impose
%                          P + margin*eD - Pop == 0 instead, eD = A'*A with
%                          A the lift of the plain term, the column of
%                          basis operators at the largest 'mult' degree of
%                          the specification (at least 1);
%            cap           double, 1 x nv, nv the registry variables of P:
%                          per direction k the bound theta_k + s_k <=
%                          cap(k) on the monomials of the lift basis, as a
%                          'subset' array on every L2 space holding all
%                          registry variables, over the degrees read or
%                          given (default []: none); the terms keep their
%                          offsets, which lower the bound where the weight
%                          of a term acts;
%          for a 'dpvar' or 'polynomial' P, OPTS is not used.
%
% OUTPUT
% - prog:  'struct', the program with the Gram matrices of Pop declared
%          positive semidefinite and P - Pop == 0 imposed; for a 'dpvar'
%          or 'polynomial' P, the program with the SOS inequality;
% - Pop:   'dopvar' object on the spaces of P, the positive operator Pop
%          ([] for a 'dpvar' or 'polynomial' P);
% - info:  'struct' with fields
%            deg         the degree specification used: a 1 x M cell of
%                        'struct' objects read from the support, or
%                        opts.deg;
%            terms       'struct' with fields codes and offsets (double
%                        rows), the terms of Pop;
%            degrees     'struct' of the degrees read from the support,
%                        with fields Dmin, Dmax, Mdeg, J, Jmax, D, w, wR,
%                        dD, dw, wmin, nv, M, vars, isR, Dmax_rule,
%                        cap_inc, or struct('given',true) with opts.deg;
%            include     1 x M cell, entry k the multi-index rows of the
%                        basis operators kept on space k ([] when not
%                        pruned);
%            direct      'struct', the construction record of Pop (field
%                        basis_list among others);
%            tag         double, the index in prog.sopeq of the record of
%                        P - Pop == 0 (0 when PROG has no field sopeq);
%            margin      double, the margin used;
%            eD          'opvar' object, eD ([] for margin 0);
%            N           double row, the Gram dimension of each term;
%            ndec        double, sum(N.^2), the entries Pop adds to
%                        prog.decvartable (Pop has N.*(N+1)/2 distinct
%                        decision variables per term);
%            sized_from  char, the source of the degrees: 'the operator',
%                        'opts.like', 'opts.deg', or, for the default
%                        like, the number of operators of prog.sop_pos
%                        found followed by ' operator(s) of poslpivar';
%            cap         double row, opts.cap ([] when none);
%          for a 'dpvar' or 'polynomial' P, struct('legacy',true).
%
% See also SOSINEQ, LPI_EQ_SOP, POSLPIVAR_DIRECT, GET_LIFT_DEGS,
%          EQ_OPTS_SOPVAR, LPIPROGRAM.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_ineq_sop
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
% Initial coding MMP, 10/08/2026

% Check the inputs; OPTS is optional
if nargin<2
    error("Not enough input arguments.")
end
if nargin<3,    opts = [];  end

% % % Not a container: 'dpvar', 'polynomial' as an SOS inequality; others
% an error
if ~(isa(P,'dopvar') || isa(P,'opvar'))
    if isa(P,'dpvar') || isa(P,'polynomial'),   prog = sosineq(prog,P);
    else
        error('lpi_ineq_sop:class','P should be an opvar, dopvar, dpvar or polynomial object; got %s.',class(P))
    end
    Pop = [];   info = struct('legacy',true);
    return
end
% OPTS as a struct; P square
if isempty(opts),   opts = struct();    end
if ~isstruct(opts)
    error("Options should be specified as a 'struct' object.")
end
if ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
    error('lpi_ineq_sop:square',"P must be square: the same spaces on both sides.")
end

% % % The degree options of the settings: a settings structure passes its
% lpi.ineq; dD, dw, wmin not given are those of prog.settings, else of the
% 'light' tier. A tier vector holds the entries for 1, 2, and 3 or more
% variables. wmin comes only with dw: a call that sets dw keeps
% w = reach + dw.
if isfield(opts,'lpi'),     tier = opts.lpi.ineq;   opts = struct();
elseif isfield(prog,'settings') && isfield(prog.settings,'lpi')
    tier = prog.settings.lpi.ineq;      % the settings of the program
else,                       tier = lpisettings_sop().lpi.ineq;  % none
end
% The term set of the tier (lpi.ineq.psatz, psatz_offset), when the tier
% has one and the call gives none: 'faces' (the 'auto' terms), 'corners'
% (the face terms and the face products), 'faces1' (the face terms one
% degree lower), 'none', 'product', or codes; a cell {1-D, 2-D, N-D} of
% these
cps = isfield(opts,'psatz');    % the call gives a term set
if isfield(tier,'psatz') && ~isfield(opts,'psatz'),    opts.psatz = tier.psatz;    end
if isfield(tier,'psatz_offset') && ~isfield(opts,'psatz_offset')
    opts.psatz_offset = tier.psatz_offset;
end
% The factor rows of the tier, only with its term set (no psatz in the
% call), so the rows never join a call's codes; a field of the call wins
for f = {'psatz_rows','psatz_rows_offset','psatz_rows_perdir'}
    if ~cps && isfield(tier,f{1}) && ~isfield(opts,f{1}),  opts.(f{1}) = tier.(f{1});  end
end
% The cross and per-space weight degrees, the lift ceiling and the
% triangular cap of the tier when the call gives none; not term options,
% so not gated on the psatz of the call
for f = {'wcross','dw_space','Dmax','cap_inc'}
    if isfield(tier,f{1}) && ~isfield(opts,f{1}),  opts.(f{1}) = tier.(f{1});  end
end
% A tier psatz as a cell {1-D, 2-D, N-D}: the entry for the registry size,
% as for the tier vectors of dD, dw, wmin below
if ~cps && isfield(opts,'psatz') && iscell(opts.psatz)
    opts.psatz = opts.psatz{min(max(numel(P.vars),1),end)};
end
% dD, dw, wmin: the entry of the tier vector for the registry size
pick = @(v) v(min(max(numel(P.vars),1),end));
if ~isfield(opts,'dD') || isempty(opts.dD),     opts.dD = pick(tier.dD);    end
if ~isfield(opts,'dw') || isempty(opts.dw)
    opts.dw = pick(tier.dw);
    if ~isfield(opts,'wmin') || isempty(opts.wmin),     opts.wmin = pick(tier.wmin);    end
end

% % % Spaces: the variable names and dimensions of each space, the domain
[sp,dm] = opvar_space_list(P,'out');
dom = struct('vars',{P.vars},'dom',P.dom);
M = numel(sp);

% % % Without opts.like and opts.deg: the support of the positive
% operators recorded in prog.sop_pos that P depends on. Read from P
% itself, a P linear in a decision operator gets a certificate for every
% value of that operator. P depends on L when it carries the decision
% names of L: up to five sampled names, each one O(q) strcmp over the
% list of P.
src = 'the operator';
hasopt = @(f) isfield(opts,f) && ~isempty(opts.(f));
if hasopt('like'),  src = 'opts.like';  end
if hasopt('deg'),   src = 'opts.deg';   end
if ~hasopt('like') && ~hasopt('deg') && isfield(prog,'sop_pos') && isa(P,'dopvar')
    zP = dvar_names(P);     like = {};
    for i = 1:numel(prog.sop_pos)
        L = prog.sop_pos{i};    zL = dvar_names(L);
        if isempty(zL) || ~all(ismember(L.vars,P.vars)),   continue,   end
        smp = zL(unique(round(linspace(1,numel(zL),min(5,numel(zL))))));
        if any(cellfun(@(z) any(strcmp(zP,z)),smp)),  like{end+1} = L;  end %#ok<AGROW>
    end
    if ~isempty(like)
        opts.like = like;
        src = sprintf('%d operator(s) of poslpivar',numel(like));
    end
end

% % % Degrees and terms: read from the support, or the given specification
if isfield(opts,'deg') && ~isempty(opts.deg)
    deg = opts.deg;     dinfo = struct('given',true);
    nv = numel(P.vars);
    % Terms: a numeric opts.psatz and its offsets, else those of the reader
    if isfield(opts,'psatz') && ~isempty(opts.psatz) && isnumeric(opts.psatz)
        codes = reshape(opts.psatz,1,[]);   offs = zeros(size(codes));
        if isfield(opts,'psatz_offset') && ~isempty(opts.psatz_offset)
            offs = reshape(opts.psatz_offset,1,[]);
            if isscalar(offs),  offs = repmat(offs,size(codes));    end
        end
        rows = zeros(0,nv);     roffs = zeros(1,0);   % rows given with the codes
        if isfield(opts,'psatz_rows') && ~isempty(opts.psatz_rows)
            rows = opts.psatz_rows;     roffs = zeros(1,size(rows,1));
            if isfield(opts,'psatz_rows_offset') && ~isempty(opts.psatz_rows_offset)
                roffs = reshape(opts.psatz_rows_offset,1,[]);
            end
        end
    else
        % The terms of the rule, read without the weight floor
        ro = opts;  ro.deg = [];
        ro.wmin = [];   % the terms of a given deg without the floor
        [~,t0,~] = get_lift_degs(P,ro);     codes = t0.codes;   offs = t0.offsets;
        rows = t0.rows;     roffs = t0.rows_offsets;
    end
    terms = struct('codes',codes,'offsets',offs,'rows',rows,'rows_offsets',roffs);
else
    % Degrees and terms read from the support
    [deg,terms,dinfo] = get_lift_degs(P,opts);
end

% % % opts.cap: per-direction bounds theta_k + s_k <= cap(k) on the lift
% basis, as a 'subset' array on every L2 space holding all registry
% variables (the only spaces that accept one), over the degrees read or
% given: subset S of [theta, s] bounded by the sum over directions of
% int_k (theta_k alone in S), mult_k (s_k alone) or min(cap_k,
% int_k+mult_k) (both). The terms keep their offsets, which lower the
% bounds where the weight of a term acts.
cap = [];
if isfield(opts,'cap') && ~isempty(opts.cap)
    cap = reshape(opts.cap,1,[]);   nv = numel(P.vars);
    if numel(cap)~=nv
        error('lpi_ineq_sop:cap',"opts.cap should have one entry per registry variable of P.")
    end
    if ~iscell(deg),    deg = repmat({deg},1,M);    end
    for k = 1:M
        % Spaces with every registry variable: int and mult per direction
        if numel(sp{k})~=nv,    continue,   end
        s = deg{k};
        if isnumeric(s) && isscalar(s),     s = struct('int',s,'mult',s);   end
        if ~isstruct(s)
            error('lpi_ineq_sop:cap',"opts.cap needs a scalar or 'struct' degree per space.")
        end
        ct = 1;     if isfield(s,'int') && ~isempty(s.int),     ct = s.int;     end
        cs = ct;    if isfield(s,'mult') && ~isempty(s.mult),   cs = s.mult;    end
        ct = reshape(ct,1,[]);  if isscalar(ct),    ct = repmat(ct,1,nv);   end
        cs = reshape(cs,1,[]);  if isscalar(cs),    cs = repmat(cs,1,nv);   end
        % The bound of each of the 2^(2nv) subsets of [theta, s]
        cj = min(cap,ct+cs);
        n = 2*nv;   sc = zeros(1,2^n);
        for i = 1:2^n
            b = bitget(i-1,1:n)>0;  bt = b(1:nv);   bs = b(nv+1:n);
            sc(i) = sum(ct.*(bt & ~bs) + cs.*(bs & ~bt) + cj.*(bt & bs));
        end
        s.subset = sc;      deg{k} = s;
    end
end

% % % The basis operators each L2 space needs: those whose diagonal term can
% reach the support of the diagonal block of P.
prune = true;
if isfield(opts,'prune') && ~isempty(opts.prune),   prune = logical(opts.prune);  end
include = cell(1,M);
if prune
    for k = 1:M
        % Spaces with variables and a stored diagonal block only
        if isempty(sp{k}) || isempty(P.C{k,k}),    continue,   end
        eo = eq_opts_sopvar(P.C{k,k});
        inc = eo.include;
        % A direction whose lower and upper cells agree is read as
        % separable, with the full-integral index 4; the lift here is not
        % separated, so 4 is expanded to the lower and upper indices 2 and
        % 3 (a superset of the unseparated reading: lossless).
        for t = 1:size(inc,2)
            r = inc(:,t)==4;
            if any(r)
                lo = inc(r,:);  up = lo;    lo(:,t) = 2;    up(:,t) = 3;
                inc = unique([inc(~r,:); lo; up],'rows');
            end
        end
        include{k} = inc;
    end
end

% % % One positive operator, all terms, then the equality.
popts = struct('psatz',terms.codes,'psatz_offset',terms.offsets);
if ~isempty(terms.rows)     % the factor rows
    popts.psatz_rows = terms.rows;  popts.psatz_rows_offset = terms.rows_offsets;
    % the per-direction offset rule of the rows (call or tier)
    if isfield(opts,'psatz_rows_perdir'),  popts.psatz_rows_perdir = opts.psatz_rows_perdir;  end
end
if prune && any(~cellfun(@isempty,include)),   popts.include = include;    end
% % % The decision variables of the certificate, counted before anything
% is declared (option 'estimate': the grids of the build below, no Gram,
% no map) and printed, so an oversized program shows before the build and
% the solve.
eopts = popts;  eopts.estimate = true;
[~,~,~,einfo] = poslpivar_direct(prog,dm,sp,dom,deg,eopts);
% Counted as prog.decvartable counts them: N^2 per Gram block of dimension
% N (the N(N+1)/2 of einfo.ndec are the distinct names of Pop).
nd = sum(einfo.N.^2);
fprintf(' - lpi_ineq_sop: %d decision variables in the certificate (Gram blocks %s, sized from %s); %d in the program.\n',...
        nd,mat2str(einfo.N),src,numel(prog.decvartable)+nd);
% Declare Pop, every term in one call
[prog,Pop,~,pinfo] = poslpivar_direct(prog,dm,sp,dom,deg,popts);
% % % The e_D margin (opts.margin = eps > 0): certify P + eps*e_D >= 0 with
% e_D = A'A the energy of the lift A of the plain term. A regularisation:
% the certified statement is the perturbed one.
eD = [];    margin = 0;
if isfield(opts,'margin') && ~isempty(opts.margin) && opts.margin>0
    margin = opts.margin;
    lopt = struct();
    if isfield(popts,'include'),    lopt.include = popts.include;   end
    dl = 1;                         % the lift degree of the plain term
    dc = deg;   if ~iscell(dc),     dc = {dc};  end
    for k = 1:numel(dc)
        if isstruct(dc{k}) && isfield(dc{k},'mult') && ~isempty(dc{k}.mult)
            dl = max([dl, reshape(dc{k}.mult,1,[])]);
        elseif isnumeric(dc{k}) && isscalar(dc{k}),  dl = max(dl,dc{k});
        end
    end
    % The lift at that degree, e_D = A_D'*A_D, and the perturbed target
    A_D = lift_opvar(dm,sp,dom,dl,lopt);
    eD = A_D'*A_D;
    P = P + margin*eD;
end
% Impose P - Pop == 0, one block per adjoint pair; TAG indexes its record
[prog,tag] = lpi_eq_sop(prog,P-Pop,'symmetric');
% The info output
info = struct('deg',{deg},'terms',terms,'degrees',dinfo,'include',{include},'direct',pinfo, ...
              'tag',tag,'margin',margin,'eD',eD,'ndec',nd,'N',einfo.N,'sized_from',src);
info.cap = cap;

end


function z = dvar_names(X)
% Z = DVAR_NAMES(X) the decision names of a 'dopvar': every 'sdopvar'
% block holds the one list of the container; {} for a fixed operator.
% INPUT  X: 'dopvar' or 'opvar' object.
% OUTPUT z: cellstr, the decision variable names ({} for an 'opvar').
z = {};
for i = 1:numel(X.C)
    if isa(X.C{i},'sdopvar'),   z = X.C{i}.Zd;  return,     end
end
end
