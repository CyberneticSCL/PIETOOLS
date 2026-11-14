function Eq = collect_eq_rows(prog,P,opts,dvars_checked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ = COLLECT_EQ_ROWS(PROG,P,OPTS,DVARS_CHECKED) returns, without adding
% them to PROG, the equations that impose P == 0 for the 'sdopvar' decision
% operator P: for each parameter cell k of P, with coefficients
% vec(C_k(d)) = A_k + B_k'*d over the decision variables d = P.Zd, the
% equations A_k(c) + B_k(:,c)'*d = 0 at the positions c in keep, those
% that the canonical multiplier form allows (the columns whose ZR monomial
% has degree 0 in every multiplier direction of cell k), as the columns of
% the block [A_k(keep).'; B_k(:,keep)]. A cell whose kept coefficients are
% all zero gives no block. With OPTS = 'symmetric', of each pair of cells
% mapped to each other by exchanging lower (2) and upper (3) integrals in
% every direction, only the one of lower linear index gives a block.
% All-zero equations are omitted. Of each set of equations that are
% multiples of one another, over all cells, only the one with the largest
% absolute coefficient is kept (the first of equals); the others are
% listed in Eq.dup. Equations without a decision variable (every
% coefficient of d at most 5e-13 |A_k(c)|) are all kept.
%
% INPUT
% - prog:   'struct', an LPI program whose prog.decvartable (a cellstr, or
%           names that cellstr(string(.)) converts) holds every name of
%           P.Zd; read only;
% - P:      'sdopvar' object in canonical multiplier form, with params.A
%           and params.B cells of 3^n3 parameters, n3 the number of
%           variables common to P.vars.in and P.vars.out;
% - opts:   (optional) 'char', 'symmetric', or [] (default []);
% - dvars_checked:  (optional) 'logical', true when every name of P.Zd is
%           known to be in prog.decvartable: the membership check is
%           skipped and Eq.pos is empty (default false).
%
% OUTPUT
% - Eq:     'struct' with fields
%             Cs     1 x nk cell, one sparse (q+1) x n_k double block per
%                    cell with equations, q = numel(P.Zd): row 1 the
%                    constant terms A_k(c), row 1+i the coefficients of
%                    Zd{i}, one column per position c of cols, in order;
%             Zd     q x 1 cellstr, the names P.Zd;
%             pos    double, q x 1, pos(i) the first row of
%                    prog.decvartable named Zd{i}; 0 x 1 when the check is
%                    skipped or Zd repeats a name;
%             cols   1 x nk cell, per block the linear positions c (double
%                    row) of its equations in the (m*NL) x (n*NR)
%                    coefficient matrix (P.dims = [m,n]; NL, NR the numbers
%                    of monomials of P.ZL, P.ZR), in ascending order;
%             cells  1 x nk double, the parameter cell (linear index into
%                    P.params.A) of each block;
%             blk    'struct' with the fields ZL, ZR, dims, vars and dom of
%                    P, NL and NR (double), psize (the size of
%                    P.params.A) and symmetric (logical, OPTS was
%                    'symmetric');
%             dup    'struct' with fields cell, col, rcell, rcol and
%                    factor, double columns with one entry per equation
%                    dropped as a multiple: its parameter cell and
%                    position, those of the equation kept for it, and f,
%                    with (dropped column) = f*(kept column) + e,
%                    |e| <= 1e-12*max|dropped column| entrywise.
%
% NOTES
% P of a class other than 'sdopvar', P not in canonical multiplier form,
% parameters of inconsistent size or count, a name of P.Zd absent from
% prog.decvartable, and OPTS other than 'symmetric' are errors.
%
% See also LPI_EQ_SDOPVAR, IMPOSE_EQ_ROWS, LPI_EQ_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - collect_eq_rows
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
% Initial coding MMP, 09/30/2026: the collect mode of 'lpi_eq_sdopvar' as
%                  its own routine (there, a second output and a 4th input
%                  selected it, so '[prog,Eq] = lpi_eq_sdopvar(...)'
%                  silently imposed nothing). The body is lpi_eq_sdopvar's
%                  from the input checks to the Eq struct, moved verbatim
%                  with its markers (MMP 08/28, 08/29, 09/26 and
%                  09/30/2026), as is the local 'canonical_positions'.
% MMP, 10/02/2026: Eq.pos, the table row of each name, from the membership
%                  check's own unique (its first-occurrence output), so
%                  'impose_eq_rows' can write prog.expr without soseq, whose
%                  getequation matched all N program names against each
%                  batch. Same verdict and errors; the N-length arrays of
%                  the check are freed before the parameter loop.
% MMP, 10/08/2026: Eq.cols, Eq.cells and Eq.blk: which coefficient position
%                  each written row constrains, per parameter, with the
%                  block's bases, so that 'lpigetdual_sop' can assemble the
%                  solver's multipliers into the dual kernel. Rows unchanged.
% MMP, 10/28/2026: Header to the header protocol, step comments in the body
%                  and the header of the local function; code unchanged.
%                  Deleted: the callers' use of the routine, a notes
%                  reference, and the reason the kept coefficients determine
%                  the operator.
% MMP, 11/01/2026: Equations that are multiples of another equation of P,
%                  over all cells, are dropped (local 'drop_multiples');
%                  Eq.dup lists them. A multiple adds no constraint; the
%                  test holds to 1e-12 of the largest coefficient, entries
%                  below 5e-13 of it (round-off of the products that build
%                  P) left out. 3-D heatNd stability LPI: 33,231 -> 23,268
%                  rows written, MOSEK 283 -> 176 s; the rows kept are
%                  unchanged bit for bit. Cs{k} holds only the columns of
%                  cols{k}. Help: description and the OUTPUT entries Cs,
%                  cols and dup.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed and
%                  in date order, (was) comment lines and a note on (was)
%                  code removed, the notes on the canonical positions, the
%                  membership check and the row layout shortened; code
%                  unchanged.
% MMP, 11/10/2026: Equations without a decision variable are all kept
%                  (drop_multiples; the 11/01 entry no longer holds for
%                  them). Any two are multiples, so the one kept was the
%                  largest |A_k(c)|: programs that differ only in A, such
%                  as the lam-affine family of heatNd_poincare (target
%                  Q0 + lam Q1), wrote different equations, and
%                  test_heatNd_poincare failed from 11/01. Other equations
%                  unchanged. Help: description and drop_multiples.

% % % Check the inputs
if isa(P,'polynomial') || isa(P,'double')
    error('Enforcing equality constraints on fixed values or polynomials is not supported.')
elseif isa(P,'sopvar')
    error("Input of type 'sopvar' carries no decision variables; use 'eq' to test whether it is zero.")
elseif ~isa(P,'sdopvar')
    error("Input must be of type 'sdopvar'.")
end
if ~isa(prog,'struct') || ~isfield(prog,'decvartable')
    error("First input must be an LPI program structure; see 'lpiprogram'.")
end

% % % Extract the structure of the operator
% dvars = cellstr(string(P.Zd(:)));                                         % MMP, 09/26/2026 (was)
% A cellstr already is the normal form; the round trip only costs O(q).     % MMP, 09/26/2026
if iscellstr(P.Zd)                                                          % MMP, 09/26/2026
    dvars = P.Zd(:);                                                        % MMP, 09/26/2026
else                                                                        % MMP, 09/26/2026
    dvars = cellstr(string(P.Zd(:)));                                       % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
% Sizes: q names; the m x n matrix; NL and NR monomials; nC coefficients
% per parameter
q = numel(dvars);
m = P.dims(1);      n = P.dims(2);
nL = cellfun(@numel,P.ZL(:)');      NL = prod([nL,1]);
nR = cellfun(@numel,P.ZR(:)');      NR = prod([nR,1]);
nC = (m*NL)*(n*NR);

% The constant (A) and decision (B) parts of the parameters, as cells
params_A = P.params.A;
params_B = P.params.B;
if ~iscell(params_A) || ~iscell(params_B)
    error("Parameters of an 'sdopvar' object should be stored as cell arrays.")
end

% The parameter cell array is indexed by the variables common to the input
% and output spaces, in sorted order, and it is those directions that can
% carry a multiplier.
vars_S3 = intersect(P.vars.in,P.vars.out);
n3 = numel(vars_S3);
if numel(params_A)~=3^n3
    error("Expected 3^n3 = "+num2str(3^n3)+" parameters for "+num2str(n3)...
          +" common variables, found "+num2str(numel(params_A))+".")
end
% Position of each shared variable among the input variables
[~,posR] = ismember(vars_S3,P.vars.in);

% Constraining only the positions the canonical form allows is equivalent to  % MMP, 08/29/2026
% constraining the operator ONLY because the rest are zero. If they are not,  % MMP, 08/29/2026
% this routine would silently impose no constraint at all where it should   % MMP, 08/29/2026
% impose one, which is a soundness failure rather than a wrong answer. The  % MMP, 08/29/2026
% constructor guarantees the invariant, so this can only fire for parameters  % MMP, 08/29/2026
% assigned directly to the properties, which bypasses it. The check inspects  % MMP, 08/29/2026
% only the forbidden positions, so it costs numel(C) per parameter and      % MMP, 08/29/2026
% nothing in the number of decision variables.                              % MMP, 08/29/2026
[is_canon,canon_info] = is_canonical_multiplier(P.params,P.vars,P.ZL,P.ZR,P.dims);% MMP, 08/29/2026
if ~is_canon                                                                % MMP, 08/29/2026
    error("The operator is not in canonical multiplier form, so constraining its "...% MMP, 08/29/2026
          +"coefficients would not constrain the operator. "+string(canon_info.message)...% MMP, 08/29/2026
          +" Rebuild it with the 'sdopvar' constructor rather than assigning to its "...% MMP, 08/29/2026
          +"properties.")                                                   % MMP, 08/29/2026
end                                                                         % MMP, 08/29/2026

% Every decision variable of the operator must be known to the program.
% The check hashes all of prog.decvartable, so a caller that has verified   % MMP, 09/26/2026
% this list already (dvars_checked) skips it.                               % MMP, 09/26/2026
% The check also gives Eq.pos, the table row of each name (see OUTPUT).     % MMP, 10/02/2026
pos = zeros(0,1);       % stays empty when the check is skipped             % MMP, 10/02/2026
% if q>0                                                                    % MMP, 09/26/2026 (was)
if q>0 && ~(nargin>=4 && dvars_checked)                                     % MMP, 09/26/2026
%   known = cellstr(string(prog.decvartable(:)));                           % MMP, 09/26/2026 (was)
    if iscellstr(prog.decvartable)                                          % MMP, 09/26/2026
        known = prog.decvartable(:);                                        % MMP, 09/26/2026
    else                                                                    % MMP, 09/26/2026
        known = cellstr(string(prog.decvartable(:)));                       % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
%   is_new = ~ismember(dvars,known);                                        % MMP, 09/26/2026 (was)
    % One unique over [known;dvars]: ic labels each distinct name, and a    % MMP, 09/26/2026
    % dvar is known iff its label also occurs among the table's (exact, as  % MMP, 09/26/2026
    % ismember, and about 2x cheaper for large lists). A single name keeps  % MMP, 09/26/2026
    % ismember, whose scalar path is one strcmp. On failure ismember runs   % MMP, 09/26/2026
    % to name the first missing variable.                                   % MMP, 09/26/2026
    if q==1                                                                 % MMP, 09/26/2026
%       is_new = ~ismember(dvars,known);                                    % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
        % The scalar path's loc is the first strcmp match, the lowest row.  % MMP, 10/02/2026
        [is_known,pos] = ismember(dvars,known);  is_new = ~is_known;        % MMP, 10/02/2026
    else                                                                    % MMP, 09/26/2026
        nk = numel(known);                                                  % MMP, 09/26/2026
%       [~,~,ic] = unique([known;dvars]);                                   % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
%       lab_known = false(max(ic),1);   lab_known(ic(1:nk)) = true;         % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
        % ia = first occurrence of each label; known comes first, so a name % MMP, 10/02/2026
        % is in the table iff ia <= nk, and ia is then its lowest row - the % MMP, 10/02/2026
        % row getequation gives it (ismember + accumarray @min).            % MMP, 10/02/2026
        [~,ia,ic] = unique([known;dvars]);                                  % MMP, 10/02/2026
        pos = ia(ic(nk+1:end));                                             % MMP, 10/02/2026
        is_new = false;                                                     % MMP, 09/26/2026
%       if ~all(lab_known(ic(nk+1:end)))                                    % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
        if any(pos>nk)                                                      % MMP, 10/02/2026
            is_new = ~ismember(dvars,known);                                % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
    if any(is_new)
        error("Decision variable '"+string(dvars{find(is_new,1)})+"' does not appear "...
              +"in the program; declare it with 'lpidecvar' before imposing constraints.")
    end
    % A name repeated in P.Zd: soseq sums its rows and drops a column that  % MMP, 10/02/2026
    % cancels, so leave such a list to soseq (empty pos). O(N), no sort.    % MMP, 10/02/2026
    if q>1                          % one name cannot repeat                % MMP, 10/02/2026
        seen = false(numel(known),1);   seen(pos) = true;                   % MMP, 10/02/2026
        if nnz(seen)<q,     pos = zeros(0,1);   end                         % MMP, 10/02/2026
    end                                                                     % MMP, 10/02/2026
    % Free the check's N-length arrays before the parameter loop; known is  % MMP, 10/02/2026
    % a new N-name cellstr when the table is not one.                       % MMP, 10/02/2026
    known = [];     ia = [];    ic = [];    seen = [];                      % MMP, 10/02/2026
end

% % % Check the options
use_symmetry = false;
if nargin>=3 && ~isempty(opts)
    if ~(ischar(opts) || isstring(opts)) || ~strcmpi(opts,'symmetric')
        error("Third input is not recognized; the only supported option is 'symmetric'.")
    end
    use_symmetry = true;
end

% sz_C = [3*ones(1,n3),1];                                                  % MMP, 09/30/2026 (was)

% Constraint blocks of the parameters, collected for one prog.expr entry.   % MMP, 10/02/2026
Cs = cell(1,numel(params_A));                                               % MMP, 09/26/2026
cols = cell(1,numel(params_A));     % positions written, per parameter      % MMP, 10/08/2026

% % % Collect the constraints, one parameter at a time                      % MMP, 09/26/2026
for k=1:numel(params_A)
    % Parameter k: an empty part as a sparse zero; skipped when all zero;
    % its size checked
    Ak = params_A{k};
    Bk = params_B{k};
    if isempty(Ak)
        Ak = sparse(nC,1);
    end
    if isempty(Bk)
        Bk = sparse(q,nC);
    end
    if nnz(Ak)==0 && nnz(Bk)==0
        continue
    end
    if numel(Ak)~=nC || size(Bk,2)~=nC || size(Bk,1)~=q
        error("Parameter "+num2str(k)+" has inconsistent dimensions.")
    end

    % Determine the multi-index of this parameter.
%   gam = ones(1,n3);                                                       % MMP, 09/30/2026 (was)
%   if n3>0                                                                 % MMP, 09/30/2026 (was)
%       idcs = cell(1,n3);                                                  % MMP, 09/30/2026 (was)
%       [idcs{:}] = ind2sub(sz_C,k);                                        % MMP, 09/30/2026 (was)
%       gam = cell2mat(idcs);                                               % MMP, 09/30/2026 (was)
%   end                                                                     % MMP, 09/30/2026 (was)
    gam = gamma_of_cell(k,n3);      % 1 x n3; 1 x 0 at n3 = 0               % MMP, 09/30/2026

    % Under self-adjointness a lower integral pairs with an upper integral.
    if use_symmetry && n3>0
        adj = gam;
        adj(gam==2) = 3;
        adj(gam==3) = 2;
%       adj_cell = num2cell(adj);                                           % MMP, 09/30/2026 (was)
%       if sub2ind(sz_C,adj_cell{:})<k                                      % MMP, 09/30/2026 (was)
        if cell_of_gamma(adj)<k                                             % MMP, 09/30/2026
            continue
        end
    end

    % Keep only the coefficients the canonical form allows to be nonzero:   % MMP, 08/29/2026
    % the rest are zero by the class invariant, so setting the kept ones to % MMP, 08/29/2026
    % zero sets the operator to zero (see 'canonicalize_multiplier').       % MMP, 08/29/2026
    keep = canonical_positions(gam,P.ZR,posR,m,n,NL,NR);                    % MMP, 08/29/2026
    A_eff = Ak(keep);                                                       % MMP, 08/29/2026
    B_eff = Bk(:,keep);                                                     % MMP, 08/29/2026
    if nnz(A_eff)==0 && nnz(B_eff)==0
        continue
    end
    n_eff = numel(keep);                                                    % MMP, 08/29/2026

    % The equations as the columns of a (q+1) x n_eff block, row 1 the
    % constant term: the layout of a 1 x n_eff variable-free 'dpvar' (one
    % block of q+1 rows per matrix ROW), which the writers take as is.
%   M = [reshape(full(A_eff),1,[]); B_eff];                                 % MMP, 09/26/2026 (was)
%   [sg,ii,vv] = find(M);                                                   % MMP, 09/26/2026 (was)
%   Cdp = sparse(sg,ii,vv,q+1,n_eff);                                       % MMP, 09/26/2026 (was)
%   Dk = dpvar(Cdp,zeros(1,0),{},dvars,[1,n_eff]);                          % MMP, 09/26/2026 (was)
%   used = find(any(B_eff,2));                                              % MMP, 09/26/2026 (was)
%   M = [reshape(full(A_eff),1,[]); B_eff(used,:)];                         % MMP, 09/26/2026 (was)
%   [sg,ii,vv] = find(M);                                                   % MMP, 09/26/2026 (was)
%   Cdp = sparse(sg,ii,vv,numel(used)+1,n_eff);                             % MMP, 09/26/2026 (was)
%   Dk = dpvar(Cdp,zeros(1,0),{},dvars(used),[1,n_eff]);                    % MMP, 09/26/2026 (was)
%   prog = soseq(prog,Dk);                                                  % MMP, 09/26/2026 (was)
    % Collect rather than impose: impose_eq_rows writes one prog.expr entry % MMP, 10/02/2026
    % per ~q nonzeros, not one per parameter.                               % MMP, 10/02/2026
    % Row 1 the constant term, row 1+i dvars{i}; sparse even if B is full.  % MMP, 09/26/2026
    Cs{k} = [sparse(reshape(A_eff,1,[])); sparse(B_eff)];                   % MMP, 09/26/2026
    % The rows the writer will produce: lpi_soseq (and soseq) drop a column % MMP, 10/08/2026
    % that is all zero, so the written rows are the kept positions whose    % MMP, 10/08/2026
    % column has a nonzero, in order (LPIGETDUAL_SOP reads them back).      % MMP, 10/08/2026
    nzc = full(any(Cs{k},1));                                               % MMP, 10/08/2026
    cols{k} = reshape(keep(nzc),1,[]);                                      % MMP, 10/08/2026
    % Columns of Cs{k} one to one with cols{k}, for drop_multiples          % MMP, 11/01/2026
    if ~all(nzc),   Cs{k} = Cs{k}(:,nzc);   end                             % MMP, 11/01/2026
end

% Equations that are multiples of another, over all cells, are dropped:     % MMP, 11/01/2026
% they constrain nothing more, and solvers pay for each row.                % MMP, 11/01/2026
[Cs,cols,dup] = drop_multiples(Cs,cols);                                    % MMP, 11/01/2026

% Skipped parameters stay [] and are dropped. Order is the loop's, which is % MMP, 09/26/2026
% the order of the former one-soseq-per-parameter expressions.              % MMP, 09/26/2026
% Eq = struct('Cs',{Cs(~cellfun(@isempty,Cs))},'Zd',{dvars});               % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
% Eq = struct('Cs',{Cs(~cellfun(@isempty,Cs))},'Zd',{dvars},'pos',{pos});   % MMP, 10/02/2026 % MMP, 10/08/2026 (was)
% cols, cells and blk describe the rows for the dual read-back: the         % MMP, 10/08/2026
% parameter index of each kept block, its written positions (linear into    % MMP, 10/08/2026
% the m*NL x n*NR coefficient matrix) and the block's bases and spaces.     % MMP, 10/08/2026
kept = ~cellfun(@isempty,Cs);                                               % MMP, 10/08/2026
blk = struct('ZL',{P.ZL},'ZR',{P.ZR},'dims',P.dims,'vars',P.vars,'dom',P.dom, ...% MMP, 10/08/2026
             'NL',NL,'NR',NR,'psize',size(params_A),'symmetric',use_symmetry);% MMP, 10/08/2026
% Eq = struct('Cs',{Cs(kept)},'Zd',{dvars},'pos',{pos}, ...                 % MMP, 10/08/2026 % MMP, 11/01/2026 (was)
%             'cols',{cols(kept)},'cells',find(kept),'blk',blk);            % MMP, 10/08/2026 % MMP, 11/01/2026 (was)
% dup: the equations dropped as multiples (see OUTPUT).                     % MMP, 11/01/2026
Eq = struct('Cs',{Cs(kept)},'Zd',{dvars},'pos',{pos}, ...                   % MMP, 10/08/2026
            'cols',{cols(kept)},'cells',find(kept),'blk',blk,'dup',dup);    % MMP, 11/01/2026

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keep = canonical_positions(gam,ZR,posR,m,n,NL,NR)
% KEEP = CANONICAL_POSITIONS(GAM,ZR,POSR,M,N,NL,NR) returns the linear
% positions in the (M*NL) x (N*NR) coefficient matrix of a parameter cell
% with multi-index GAM that the canonical multiplier form allows: those
% whose column has a right monomial of degree 0 in every direction GAM
% marks as a multiplier (gam(t) = 1).
% INPUT  gam:    1 x n3 double, entries 1 (multiplier), 2, 3 (integrals);
%        ZR:     cell of double columns, the exponents of the right
%                monomials in each input variable;
%        posR:   1 x n3 double, the index in ZR of shared variable t;
%        m, n:   double, the matrix dimensions of the operator;
%        NL, NR: double, the numbers of left and right monomials.
% OUTPUT keep:   double column, the positions in ascending order.

% Columns allowed: in each multiplier direction the right monomials of
% degree 0 in that variable, every monomial in the other variables
nR = cellfun(@numel,ZR(:)');
okB = true(NR,1);
for t = find(gam==1)
    p = posR(t);
    v = (ZR{p}(:)==0);
    kk = true(1,1);
    for i = 1:numel(nR)
        if i==p
            kk = kron(kk,v);
        else
            kk = kron(kk,true(nR(i),1));
        end
    end
    okB = okB & kk;
end

% Monomials are ordered as ZR(x) = x_1^ZR{1} o ... o x_N^ZR{N}, so a column
% of the coefficient matrix is (matrix column, monomial) with the monomial
% as the inner index.
nrow = m*NL;
[bG,jG,rG] = ndgrid(find(okB)-1,(0:n-1)',(0:nrow-1)');
keep = sort((jG(:)*NR + bG(:))*nrow + rG(:) + 1);

end



% BEGIN MMP, 11/01/2026: local function drop_multiples, new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cs,cols,dup] = drop_multiples(Cs,cols)
% [CS,COLS,DUP] = DROP_MULTIPLES(CS,COLS) removes from the blocks CS every
% column that is a multiple of another column of CS. Entries of at most
% 5e-13 times the largest absolute entry of their column are left out of
% the comparison. Each column is scaled by s, the sign of its first entry
% compared times its largest absolute entry; columns with the same rows of
% entries compared, and scaled entries there that agree to 1e-12, form a
% group, and of each group only the column of largest |s| (the first of
% equals) is kept. COLS is pruned with CS. Columns whose only entry
% compared is in row 1, the constant term, are all kept.
% INPUT  Cs:   1 x nk cell; Cs{k} [] or a sparse (q+1) x L_k 'double'
%              block with no all-zero column;
%        cols: 1 x nk cell; cols{k} the L_k positions of the columns of
%              Cs{k}, a 'double' row.
% OUTPUT Cs, cols: the inputs without the dropped columns;
%        dup:  'struct' with fields cell, col, rcell, rcol and factor,
%              'double' columns, one entry per dropped column: its block k
%              and position, those of the column kept for it, and f with
%              (dropped column) = f*(kept column) + e, |e| <= 1e-12 times
%              the largest absolute entry of the dropped column.

tol = 1e-12;
dup = struct('cell',zeros(0,1),'col',zeros(0,1),'rcell',zeros(0,1), ...
             'rcol',zeros(0,1),'factor',zeros(0,1));
% Blocks with columns; the columns of block ks(t) are o(t)+1:o(t+1) in the
% joint numbering
ks = find(~cellfun(@isempty,Cs));
L = cellfun(@(C) size(C,2),Cs(ks));
o = cumsum([0, L]);     nt = o(end);
if nt<2,    return,     end

% Pass 1, per column: the entries compared, those above thr = tol/2 times
% the largest absolute entry (smaller ones are round-off of cancelled sums
% in the products that build P, which change the pattern of a multiple);
% their number, first row, the scale s and two weighted sums of the scaled
% entries (the hash). The weights are Weyl sequences of the row index,
% evaluated at the nonzeros only: no array of length q. find lists a block
% by column, rows ascending, so fp is the first entry of each column.
% ~(<=) keeps a NaN entry, so that its column is never matched.
cnt = zeros(nt,1);  fr = zeros(nt,1);   s = zeros(nt,1);    h = zeros(nt,2);
thr = zeros(nt,1);
for t = 1:numel(ks)
    [ii,jj,vv] = find(Cs{ks(t)});
    ii = ii(:);     jj = jj(:);     vv = vv(:);     % rows for a 1 x L block (q = 0)
    J = o(t)+(1:L(t))';
    mx = accumarray(jj,abs(vv),[L(t),1],@max);
    thr(J) = (tol/2)*mx;
    sg = ~(abs(vv) <= thr(J(jj)));
    ii = ii(sg);    jj = jj(sg);    vv = vv(sg);
    c = accumarray(jj,1,[L(t),1]);
    fp = cumsum([1; c(1:end-1)]);
    st = sign(vv(fp)).*mx;
    vv = vv./st(jj);
    cnt(J) = c;     fr(J) = ii(fp);     s(J) = st;
    h(J,1) = accumarray(jj,vv.*(mod(ii*0.6180339887498949,1)-0.5),[L(t),1]);
    h(J,2) = accumarray(jj,vv.*(mod(ii*0.4142135623730950,1)-0.5),[L(t),1]);
end

% Rounds over the columns act still to be grouped, at first all of them.
% Candidate groups: equal count, first row and hashes rounded to 1e-8, as
% multiples have equal scaled columns. The representative of a group is
% its member of largest |s|, so that a solver that zeroes entries below an
% absolute tolerance (MOSEK: 1e-12) cannot erase the kept row while a
% larger multiple of it is dropped. A hash match is no proof: every
% candidate is compared with its representative entrywise. A near multiple
% can share the key and, as representative, refuse the true multiples of
% its group, so the refused candidates are grouped again among themselves.
key = [cnt, fr, round(h*1e8)];
act = (1:nt)';      drop = zeros(0,1);      r = zeros(0,1);     gathered = false;
% Columns whose only entry compared is row 1 (no decision variable) stay:   % MMP, 11/10/2026
% any two are multiples, so the one kept would follow the constant terms,   % MMP, 11/10/2026
% and programs that differ only in those would write different equations.   % MMP, 11/10/2026
act = act(~(cnt==1 & fr==1));                                               % MMP, 11/10/2026
while numel(act)>=2
    % Groups of act; cand the members that are not representatives, rc
    % their representatives
    [~,~,g] = unique(key(act,:),'rows');
    [~,ord] = sortrows([g, -abs(s(act)), act]);
    first = [true; diff(g(ord))~=0];
    rg = act(ord(first));
    rp = zeros(numel(act),1);   rp(ord) = rg(cumsum(first));
    cand = act(rp~=act);    rc = rp(rp~=act);
    if isempty(cand),   break,  end
    % Pass 2, first round only (later candidates and representatives are
    % among the first candidates): rows R and scaled values V of the
    % candidates and representatives, column after column, column j from
    % ptr(j) on
    if ~gathered
        use = false(nt,1);  use(cand) = true;   use(rc) = true;
        ci = cnt.*use;
        ptr = cumsum([1; ci(1:end-1)]);
        R = zeros(sum(ci),1);   V = zeros(sum(ci),1);
        for t = 1:numel(ks)
            J = o(t)+(1:L(t))';
            j0 = J(find(use(J),1));
            if isempty(j0),     continue,   end
            [ii,jj,vv] = find(Cs{ks(t)});
            ii = ii(:);     jj = jj(:);     vv = vv(:);
            sel = use(J(jj)) & ~(abs(vv) <= thr(J(jj)));     % as in pass 1
            p = ptr(j0) + (0:nnz(sel)-1)';
            R(p) = ii(sel);     V(p) = vv(sel)./s(J(jj(sel)));
        end
        gathered = true;
    end
    % Entry l of a candidate against entry l of its representative (the
    % key makes the counts equal): the same row, scaled values within tol;
    % ~(<=) so that a NaN fails. grp: the candidate of each compared entry,
    % a column (repelem returns a row for one candidate).
    nc = cnt(cand);
    grp = repelem((1:numel(cand))',nc);     grp = grp(:);
    st0 = cumsum([0; nc]);
    off = (1:st0(end))' - st0(grp) - 1;
    ia = ptr(cand(grp)) + off;
    ib = ptr(rc(grp)) + off;
    bad = R(ia)~=R(ib) | ~(abs(V(ia)-V(ib))<=tol);
    ok = accumarray(grp,double(bad),[numel(cand),1])==0;
    % Few rounds (one when no candidate is refused), O(nt) each
    drop = [drop; cand(ok)];    r = [r; rc(ok)];  %#ok<AGROW>
    act = cand(~ok);
end
if isempty(drop),   return,     end

% The dropped columns with the columns kept for them, then the blocks and
% positions without the dropped columns
tb = repelem((1:numel(ks))',L(:));      % block (index into ks) of each column
tb = tb(:);
pos = [cols{ks}].';
dup = struct('cell',reshape(ks(tb(drop)),[],1),'col',pos(drop), ...
             'rcell',reshape(ks(tb(r)),[],1),'rcol',pos(r),'factor',s(drop)./s(r));
keep = true(nt,1);  keep(drop) = false;
for t = 1:numel(ks)
    kt = keep(o(t)+1:o(t+1));
    if ~all(kt)
        Cs{ks(t)} = Cs{ks(t)}(:,kt);
        cols{ks(t)} = cols{ks(t)}(kt);
    end
end

end
% END MMP, 11/01/2026
