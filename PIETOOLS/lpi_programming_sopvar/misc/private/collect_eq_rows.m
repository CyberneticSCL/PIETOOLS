function Eq = collect_eq_rows(prog,P,opts,dvars_checked)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EQ = COLLECT_EQ_ROWS(PROG,P,OPTS,DVARS_CHECKED) returns, without imposing
% them, the equality constraints that enforce P==0 for an 'sdopvar' P: one
% row block per constrained parameter cell. 'lpi_eq_sdopvar' is this
% followed by 'impose_eq_rows'; 'lpi_eq_cdopvar' collects the blocks of
% every container block over one decision list and imposes them together.
%
% The rows are the stored coefficients of each parameter, restricted to
% the positions the canonical multiplier form allows to be nonzero: by the
% class invariant the operator vanishes exactly when they do (see
% 'lpi_eq_sdopvar'; sopvar_implementation_notes.pdf Sec. 8.1.6 for
% equality by coefficients).
%
% INPUT
% - prog:   'struct' LPI program; every entry of P.Zd must be in
%           prog.decvartable. Read only: prog is not modified;
% - P:      'sdopvar' object;
% - opts:   [] or 'symmetric', as for 'lpi_eq_sdopvar';
% - dvars_checked: true if the caller has already verified that every entry
%           of P.Zd is in prog.decvartable, so the check, which hashes all
%           of prog.decvartable, is skipped. 'lpi_eq_cdopvar' uses it to
%           check once per distinct list rather than once per block;
%
% OUTPUT
% - Eq:     struct with fields Cs, one sparse (q+1) x n_k block per
%           constrained parameter (row 1 the constant term, row 1+i
%           decision variable Eq.Zd{i}), in the order one soseq per
%           parameter imposed them, and Zd, the q x 1 cellstr of names.
%           Field pos (q x 1): pos(i) = first row of prog.decvartable named % MMP, 10/02/2026
%           Zd{i}, from the membership check; 0 x 1 when the check is       % MMP, 10/02/2026
%           skipped (dvars_checked) or Zd repeats a name. 'impose_eq_rows'  % MMP, 10/02/2026
%           writes the rows directly ('lpi_soseq') when it has pos.         % MMP, 10/02/2026
%
% See also LPI_EQ_SDOPVAR, IMPOSE_EQ_ROWS, LPI_EQ_CDOPVAR.

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
%                  its own routine. 'lpi_eq_sdopvar' selected collect-only
%                  by a second output and a 4th input, and impose-only by a
%                  struct P, so '[prog,Eq] = lpi_eq_sdopvar(...)' silently
%                  imposed nothing. The body below is lpi_eq_sdopvar's from
%                  the input checks to the Eq struct, moved verbatim with
%                  its markers (MMP 08/28, 08/29, 09/26 and 09/30/2026), as
%                  is the local 'canonical_positions'. Not moved: the
%                  nargout switch that imposed the rows ('lpi_eq_sdopvar'
%                  calls 'impose_eq_rows' instead). One comment renames
%                  impose_rows -> impose_eq_rows (marked below).
% MMP, 10/08/2026: Eq.cols, Eq.cells and Eq.blk: which coefficient position
%                  each written row constrains, per parameter, with the
%                  block's bases, so that 'lpigetdual_sop' can assemble the
%                  solver's multipliers into the dual kernel. Rows unchanged.
% MMP, 10/02/2026: Eq.pos, the table row of each name, from the membership
%                  check's own unique (its first-occurrence output), so
%                  'impose_eq_rows' can write prog.expr without soseq,
%                  whose getequation matched all N program names against
%                  each batch (0.54 s per call at N = 2.7e6, measured).
%                  Same verdict and errors; cost +7 B x N transient (ia
%                  replaces a logical), 8 B x q held, and an O(N) repeat
%                  check (1 B x N, q > 1 only). The check's N-length arrays
%                  (known, ia, ic) are now freed before the parameter loop;
%                  they used to live through it.

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
q = numel(dvars);
m = P.dims(1);      n = P.dims(2);
nL = cellfun(@numel,P.ZL(:)');      NL = prod([nL,1]);
nR = cellfun(@numel,P.ZR(:)');      NR = prod([nR,1]);
nC = (m*NL)*(n*NR);

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
[~,posR] = ismember(vars_S3,P.vars.in);

% Constraining only the positions the canonical form allows is equivalent to  % MMP, 08/29/2026
% constraining the operator ONLY because the rest are zero. If they are not,  % MMP, 08/29/2026
% this routine would silently impose no constraint at all where it should     % MMP, 08/29/2026
% impose one, which is a soundness failure rather than a wrong answer. The    % MMP, 08/29/2026
% constructor guarantees the invariant, so this can only fire for parameters  % MMP, 08/29/2026
% assigned directly to the properties, which bypasses it. The check inspects  % MMP, 08/29/2026
% only the forbidden positions, so it costs numel(C) per parameter and        % MMP, 08/29/2026
% nothing in the number of decision variables.                                % MMP, 08/29/2026
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
    % dvar is known iff its label also occurs among the table's. Exact:     % MMP, 09/26/2026
    % unique and ismember both match names char for char. Measured ~2x      % MMP, 09/26/2026
    % cheaper than ismember for 2 to 3.8e5 names against 7.6e5; a single    % MMP, 09/26/2026
    % name keeps ismember, whose scalar path is one strcmp (18 vs 171 ms).  % MMP, 09/26/2026
    % On failure ismember runs to name the first missing variable.          % MMP, 09/26/2026
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

% Constraint blocks of the parameters, collected for one soseq (see Eq).    % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
% Constraint blocks of the parameters, collected for one prog.expr entry.   % MMP, 10/02/2026
Cs = cell(1,numel(params_A));                                               % MMP, 09/26/2026
cols = cell(1,numel(params_A));     % positions written, per parameter      % MMP, 10/08/2026

% % % Impose the constraints, one parameter at a time                       % MMP, 09/26/2026 (was)
% % % Collect the constraints, one parameter at a time                      % MMP, 09/26/2026
for k=1:numel(params_A)
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

    % Keep only the coefficients the canonical form allows to be nonzero. In  % MMP, 08/29/2026
    % a multiplier direction the rest are zero by the class invariant, so     % MMP, 08/29/2026
    % constraining them would add nothing; and because the invariant holds,   % MMP, 08/29/2026
    % setting the kept coefficients to zero is equivalent to setting the      % MMP, 08/29/2026
    % operator to zero. Earlier versions had to contract the multiplier       % MMP, 08/29/2026
    % directions here instead, since the coefficients were not then           % MMP, 08/29/2026
    % determined by the operator; see 'canonicalize_multiplier'.              % MMP, 08/29/2026
    keep = canonical_positions(gam,P.ZR,posR,m,n,NL,NR);                     % MMP, 08/29/2026
    A_eff = Ak(keep);                                                        % MMP, 08/29/2026
    B_eff = Bk(:,keep);                                                      % MMP, 08/29/2026
    if nnz(A_eff)==0 && nnz(B_eff)==0
        continue
    end
    n_eff = numel(keep);                                                     % MMP, 08/29/2026

    % Wrap the affine expressions as a variable-free 'dpvar' and impose the
    % equality. The constraints are laid out as a ROW of n_eff scalar
    % entries rather than a column: a dpvar coefficient matrix carries one
    % block of q+1 rows per matrix ROW, so a column layout would give it
    % n_eff*(q+1) rows, and the generic combine/compress pass inside
    % sosconstr scales with that row count. As a row the coefficient matrix
    % is (q+1) by n_eff with the same nonzeros, which is orders of magnitude
    % cheaper to process and yields byte-identical At, b and Z.
%   M = [reshape(full(A_eff),1,[]); B_eff];                                 % MMP, 09/26/2026 (was)
%   [sg,ii,vv] = find(M);                                                   % MMP, 09/26/2026 (was)
%   Cdp = sparse(sg,ii,vv,q+1,n_eff);                                       % MMP, 09/26/2026 (was)
%   Dk = dpvar(Cdp,zeros(1,0),{},dvars,[1,n_eff]);                          % MMP, 09/26/2026 (was)
    % Pass only the decision variables this parameter uses. sosconstr's     % MMP, 09/26/2026
    % combine/compress drop the rest anyway, but only after O(q) passes     % MMP, 09/26/2026
    % over all q names (char + sortrows in DPVuniquedvar). In the 2-D Hinf  % MMP, 09/26/2026
    % build a parameter uses 1 to 1.0e5 of q = 3.8e5. At, b, Z unchanged:   % MMP, 09/26/2026
    % compress would remove exactly these zero rows, and getequation places % MMP, 09/26/2026
    % the remaining rows by name.                                           % MMP, 09/26/2026
%   used = find(any(B_eff,2));                                              % MMP, 09/26/2026 (was)
%   M = [reshape(full(A_eff),1,[]); B_eff(used,:)];                         % MMP, 09/26/2026 (was)
%   [sg,ii,vv] = find(M);                                                   % MMP, 09/26/2026 (was)
%   Cdp = sparse(sg,ii,vv,numel(used)+1,n_eff);                             % MMP, 09/26/2026 (was)
%   Dk = dpvar(Cdp,zeros(1,0),{},dvars(used),[1,n_eff]);                    % MMP, 09/26/2026 (was)
%   prog = soseq(prog,Dk);                                                  % MMP, 09/26/2026 (was)
    % Collect rather than impose: joined soseq calls (impose_eq_rows) pay   % MMP, 09/30/2026 % MMP, 10/02/2026 (was)
    % getequation's scan of the q program names per ~q nonzeros, not per    % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
    % parameter.                                                            % MMP, 09/26/2026 % MMP, 10/02/2026 (was)
    % Collect rather than impose: impose_eq_rows writes one prog.expr entry % MMP, 10/02/2026
    % per ~q nonzeros, not one per parameter.                               % MMP, 10/02/2026
    % Row 1 the constant term, row 1+i dvars{i}; sparse even if B is full.  % MMP, 09/26/2026
    Cs{k} = [sparse(reshape(A_eff,1,[])); sparse(B_eff)];                   % MMP, 09/26/2026
    % The rows the writer will produce: lpi_soseq (and soseq) drop a column % MMP, 10/08/2026
    % that is all zero, so the written rows are the kept positions whose    % MMP, 10/08/2026
    % column has a nonzero, in order (LPIGETDUAL_SOP reads them back).      % MMP, 10/08/2026
    nzc = full(any(Cs{k},1));                                               % MMP, 10/08/2026
    cols{k} = reshape(keep(nzc),1,[]);                                      % MMP, 10/08/2026
end

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
Eq = struct('Cs',{Cs(kept)},'Zd',{dvars},'pos',{pos}, ...                   % MMP, 10/08/2026
            'cols',{cols(kept)},'cells',find(kept),'blk',blk);              % MMP, 10/08/2026

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keep = canonical_positions(gam,ZR,posR,m,n,NL,NR)
% KEEP = CANONICAL_POSITIONS(...) returns the linear positions within vec(C)
% that the canonical multiplier form allows a parameter with multi-index GAM
% to use: those whose right monomial has degree 0 in every direction GAM
% marks as a multiplier. Every degree is kept in the remaining directions.
%
% The parameters of an 'sdopvar' are canonical by class invariant, so the
% positions this drops are zero and constraining them would add nothing.

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
