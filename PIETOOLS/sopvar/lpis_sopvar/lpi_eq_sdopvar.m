function [prog,Eq] = lpi_eq_sdopvar(prog,P,opts,dvars_checked)              % MMP, 09/26/2026
% function prog = lpi_eq_sdopvar(prog,P,opts,dvars_checked)                 % MMP, 09/26/2026 (was)
% function prog = lpi_eq_sdopvar(prog,P,opts)                               % MMP, 09/26/2026 (was)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROG = LPI_EQ_SDOPVAR(PROG,P) takes an LPI optimization program structure
% 'prog' and an 'sdopvar' decision variable P, and adds equality constraints
% enforcing P==0. It is the 'sdopvar' counterpart of 'lpi_eq', and applies in
% particular to the operators returned by 'possopvar'.
%
% An 'sdopvar' object represents the operator
%
%   (P x)(s) = sum_gamma int_{s'} K_gamma(s,s') I_gamma(s-s') x(s') ds'
%
%   K_gamma(s,s') = (I_m o ZL(s))^T unvec(A_gamma + B_gamma'*d) (I_n o ZR(s'))
%
% for gamma in {1,2,3}^n3, where index 1 denotes a multiplier (a factor
% delta(s_k-s_k')), 2 a lower integral and 3 an upper integral.
%
% Since the 3^n3 indicator patterns I_gamma act on disjoint regions, the
% operator vanishes exactly when every term vanishes, and a term vanishes
% exactly when its coefficients do. The second half of that relies on the
% canonical multiplier form, which every 'sdopvar' satisfies: in a direction
% where gamma_k=1 the factor delta(s_k-s_k') identifies s_k with s_k', so
% degrees could otherwise be moved freely between ZL and ZR without changing
% the operator, and constraining the stored coefficients would be strictly
% stronger than constraining the operator. The canonical form pins that
% freedom by keeping all of the degree on the left, so the coefficients are
% determined by the operator and constraining them is exactly right. See
% 'canonicalize_multiplier' and the class documentation.
%
% This routine therefore constrains the stored coefficients directly, after
% dropping the positions the canonical form requires to be zero.
%
% INPUT
% - prog:   'struct' specifying an LPI program structure to modify (see also
%           'lpiprogram'). All decision variables of P must already appear in
%           'prog.decvartable';
% - P:      'sdopvar' object of which to enforce P==0. Parameters that are
%           empty or identically zero generate no constraints;
%           (internal) or a batch EQ returned by earlier calls, which is    % MMP, 09/26/2026
%           then imposed as it stands; see EQ below;                        % MMP, 09/26/2026
% - opts:   (optional) set opts = 'symmetric' if P is known to be
%           self-adjoint, to constrain only one parameter per adjoint pair.
%           Taking the adjoint exchanges lower and upper integrals in every
%           direction, so the parameters pair up under the index map
%           1->1, 2->3, 3->2;
% - dvars_checked: (optional, internal) true if the caller has already      % MMP, 09/26/2026
%           verified that every entry of P.Zd is in prog.decvartable, so    % MMP, 09/26/2026
%           the check here, which hashes all of prog.decvartable, is        % MMP, 09/26/2026
%           skipped. 'lpi_eq_cdopvar' uses it to check once per distinct    % MMP, 09/26/2026
%           list rather than once per block. Default false;                 % MMP, 09/26/2026
%
% OUTPUT
% - prog:   Same LPI program structure as the input, but with the added
%           constraints enforcing P==0.
% - Eq:     (optional, internal) if requested, the constraints are NOT      % MMP, 09/26/2026
%           imposed; prog is returned unchanged and Eq holds them: Eq.Cs,   % MMP, 09/26/2026
%           one sparse (q+1) x n_k block per constrained parameter (row 1   % MMP, 09/26/2026
%           the constant term, row 1+i decision variable Eq.Zd{i}), in the  % MMP, 09/26/2026
%           order they would be imposed. 'lpi_eq_cdopvar' joins the blocks  % MMP, 09/26/2026
%           of all its blocks and imposes them with one call;               % MMP, 09/26/2026
%
% NOTES
% To enforce P==Q, call LPI_EQ_SDOPVAR(PROG,P-Q); '@sdopvar/plus' aligns the
% monomial and decision variable bases of the two operators.
%
% The constraints of one call are joined into as few 'soseq' calls as a cap % MMP, 09/26/2026
% of numel(prog.decvartable) nonzeros each allows, so prog.expr gets a few  % MMP, 09/26/2026
% entries per call, not one per parameter. The rows sossolve assembles from % MMP, 09/26/2026
% prog.expr (At and b, in order) are the same either way.                   % MMP, 09/26/2026
%
% See also LPI_EQ, LPI_EQ_NDOPVAR, POSSOPVAR, SOSEQ, SDOPVAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - lpi_eq_sdopvar
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
% MMP, 08/28/2026: Initial coding
% MMP, 08/29/2026: Drop the multiplier contraction. The canonical multiplier
%                  form is now a class invariant, enforced by the 'sdopvar'
%                  constructor, so the stored coefficients are determined by
%                  the operator and can be constrained directly. What remains
%                  is a selection of the positions the form allows to be
%                  nonzero, which is pure indexing.
% MMP, 09/26/2026: Name handling in O(q) per call removed from the hot path.
%                  (1) Use P.Zd and prog.decvartable as they are when they are
%                  already cellstr, instead of a cellstr(string(.)) round trip
%                  over q names each; (2) optional 4th input lets
%                  'lpi_eq_cdopvar' skip the membership check for blocks whose
%                  list it has already verified; (3) hand soseq only the
%                  decision variables a parameter uses, so sosconstr's
%                  combine/compress no longer sort all q names per call.
%                  Measured, 2-D container Hinf build (q = 7.6e5, 6 calls,
%                  14 soseq): the check and conversions took 3.7 s and
%                  combine/compress 1.6 s; same program, bit for bit.
% MMP, 09/26/2026: Membership check by one unique over [table; names]
%                  instead of ismember(names,table). The reverse scan
%                  getequation took today (ismember(table,names)) does not
%                  help here: the names are half the table (3.8e5 of 7.6e5),
%                  so both sides are large; measured 261 vs 275 ms. The
%                  unique route, same process vs HEAD: 2-D Hinf 411 -> 212
%                  ms, q = 1e6 over 3 variables 461 -> 191 ms; peak +8 MB
%                  at 1.1e6 names (6.5 -> 14.5 MB). A single name keeps
%                  ismember (a strcmp there). Same verdict, same error.
% MMP, 09/26/2026: One soseq per ~q nonzeros instead of one per parameter.
%                  Every soseq scans all q names of prog.decvartable in
%                  getequation (ismember(decvartable,dvarname)), and each
%                  call hashes its names again. Measured, 2-D container Hinf
%                  (q = 7.6e5): 14 calls, 28 ms fixed each (1-entry soseq),
%                  738k names hashed for a 3.8e5 list; soseq 1.74 of the
%                  2.29 s of S6b. The parameters' columns are now collected
%                  in the old call order and imposed as row dpvars joined up
%                  to q nonzeros each (impose_rows). A second output returns
%                  them unimposed and P = that output imposes them, so
%                  'lpi_eq_cdopvar' joins all its blocks. The cap, not one
%                  soseq for everything: soseq holds ~9 copies of its input,
%                  so one uncapped call peaked at 3.4x the per-parameter
%                  memory. Measured 2-D Hinf, q = 7.6e5 / 2.3e6, old -> new:
%                  2.2 -> 1.1 s / 9.5 -> 5.0 s; peak 90 -> 138 MB / 291 ->
%                  447 MB (stock lpi_eq_2d 1.5 s, 174 MB / 6.3 s, 540 MB).
%                  Item (3) of the first 09/26 entry now applies per soseq:
%                  the names used by any of its parameters. prog.expr has
%                  fewer, longer entries; the At, b, c, K sossolve builds
%                  are bit-identical, row order included. Cost: soseq calls
%                  ~ nnz/q + name-list changes, instead of the parameter
%                  count (3^n3 per block, multiplicative in dimension);
%                  memory ~9 x 16 B x min(q + largest parameter, nnz) per
%                  soseq plus one 16 B/nnz copy of the collected blocks.
%                  Measured in 1-D and 2-D only.


% Internal batch form: P is the Eq output of earlier calls, possibly        % MMP, 09/26/2026
% several joined over one name list by 'lpi_eq_cdopvar'. Its names were     % MMP, 09/26/2026
% verified when it was collected.                                           % MMP, 09/26/2026
if isstruct(P) && isfield(P,'Cs') && isfield(P,'Zd')                        % MMP, 09/26/2026
    prog = impose_rows(prog,P);                                             % MMP, 09/26/2026
    return                                                                  % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026

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
        is_new = ~ismember(dvars,known);                                    % MMP, 09/26/2026
    else                                                                    % MMP, 09/26/2026
        nk = numel(known);                                                  % MMP, 09/26/2026
        [~,~,ic] = unique([known;dvars]);                                   % MMP, 09/26/2026
        lab_known = false(max(ic),1);   lab_known(ic(1:nk)) = true;         % MMP, 09/26/2026
        is_new = false;                                                     % MMP, 09/26/2026
        if ~all(lab_known(ic(nk+1:end)))                                    % MMP, 09/26/2026
            is_new = ~ismember(dvars,known);                                % MMP, 09/26/2026
        end                                                                 % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
    if any(is_new)
        error("Decision variable '"+string(dvars{find(is_new,1)})+"' does not appear "...
              +"in the program; declare it with 'lpidecvar' before imposing constraints.")
    end
end

% % % Check the options
use_symmetry = false;
if nargin>=3 && ~isempty(opts)
    if ~(ischar(opts) || isstring(opts)) || ~strcmpi(opts,'symmetric')
        error("Third input is not recognized; the only supported option is 'symmetric'.")
    end
    use_symmetry = true;
end

sz_C = [3*ones(1,n3),1];

% Constraint blocks of the parameters, collected for one soseq (see Eq).    % MMP, 09/26/2026
Cs = cell(1,numel(params_A));                                               % MMP, 09/26/2026

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
    gam = ones(1,n3);
    if n3>0
        idcs = cell(1,n3);
        [idcs{:}] = ind2sub(sz_C,k);
        gam = cell2mat(idcs);
    end

    % Under self-adjointness a lower integral pairs with an upper integral.
    if use_symmetry && n3>0
        adj = gam;
        adj(gam==2) = 3;
        adj(gam==3) = 2;
        adj_cell = num2cell(adj);
        if sub2ind(sz_C,adj_cell{:})<k
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
    % Collect rather than impose: joined soseq calls (impose_rows) pay      % MMP, 09/26/2026
    % getequation's scan of the q program names per ~q nonzeros, not per    % MMP, 09/26/2026
    % parameter.                                                            % MMP, 09/26/2026
    % Row 1 the constant term, row 1+i dvars{i}; sparse even if B is full.  % MMP, 09/26/2026
    Cs{k} = [sparse(reshape(A_eff,1,[])); sparse(B_eff)];                   % MMP, 09/26/2026
end

% Skipped parameters stay [] and are dropped. Order is the loop's, which is % MMP, 09/26/2026
% the order of the former one-soseq-per-parameter expressions.              % MMP, 09/26/2026
Eq = struct('Cs',{Cs(~cellfun(@isempty,Cs))},'Zd',{dvars});                 % MMP, 09/26/2026
if nargout<2                                                                % MMP, 09/26/2026
    prog = impose_rows(prog,Eq);                                            % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % MMP, 09/26/2026
function prog = impose_rows(prog,Eq)                                        % MMP, 09/26/2026
% PROG = IMPOSE_ROWS(PROG,EQ) imposes the blocks EQ.Cs{:}, each (q+1) x n_k % MMP, 09/26/2026
% over the names EQ.Zd with row 1 the constant term, joining consecutive    % MMP, 09/26/2026
% blocks into as few soseq calls as the cap below allows, each on one       % MMP, 09/26/2026
% variable-free row dpvar (see the layout note in the main loop). Columns   % MMP, 09/26/2026
% keep their order, so the At and b produced are those of one soseq per     % MMP, 09/26/2026
% block, concatenated as sossolve concatenates prog.expr.                   % MMP, 09/26/2026
% Cap: at most N = numel(prog.decvartable) nonzeros per soseq; a larger     % MMP, 09/26/2026
% block goes alone. Each soseq pays an O(N) scan of the program's names     % MMP, 09/26/2026
% (getequation), which N nonzeros of work amortize; and soseq holds ~9      % MMP, 09/26/2026
% copies of its coefficients (combine, compress, getequation), so each      % MMP, 09/26/2026
% soseq's transient is O(N) bytes. Eq.Cs itself, one copy of all the        % MMP, 09/26/2026
% collected coefficients (O(nnz)), lives until the last soseq. Measured,    % MMP, 09/26/2026
% 2-D container Hinf, q = 2.3e6: one uncapped soseq peaked at 974 MB,       % MMP, 09/26/2026
% capped 447 MB (per parameter 291 MB, stock lpi_eq_2d 540 MB).             % MMP, 09/26/2026
if isempty(Eq.Cs),  return,     end                                         % MMP, 09/26/2026
cap = numel(prog.decvartable);                                              % MMP, 09/26/2026
nz = cellfun(@nnz,Eq.Cs);                                                   % MMP, 09/26/2026
k0 = 1;                                                                     % MMP, 09/26/2026
while k0<=numel(Eq.Cs)                                                      % MMP, 09/26/2026
    k1 = k0;    tot = nz(k0);                                               % MMP, 09/26/2026
    while k1<numel(Eq.Cs) && tot+nz(k1+1)<=cap                              % MMP, 09/26/2026
        k1 = k1+1;  tot = tot+nz(k1);                                       % MMP, 09/26/2026
    end                                                                     % MMP, 09/26/2026
    M = [Eq.Cs{k0:k1}];                                                     % MMP, 09/26/2026
    % Keep the constant row and the rows of names in use; compress would    % MMP, 09/26/2026
    % drop the others only after O(q) passes over their names.              % MMP, 09/26/2026
    rows = find(any(M,2));                                                  % MMP, 09/26/2026
    rows = [1; rows(rows>1)];                                               % MMP, 09/26/2026
    Dk = dpvar(M(rows,:),zeros(1,0),{},Eq.Zd(rows(2:end)-1),[1,size(M,2)]); % MMP, 09/26/2026
    % Free the joined copy before soseq makes its own.                      % MMP, 09/26/2026
    M = [];                                                                 % MMP, 09/26/2026
    prog = soseq(prog,Dk);                                                  % MMP, 09/26/2026
    k0 = k1+1;                                                              % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026
end                                                                         % MMP, 09/26/2026



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
