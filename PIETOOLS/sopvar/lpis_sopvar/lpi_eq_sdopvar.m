function prog = lpi_eq_sdopvar(prog,P,opts)
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
% - opts:   (optional) set opts = 'symmetric' if P is known to be
%           self-adjoint, to constrain only one parameter per adjoint pair.
%           Taking the adjoint exchanges lower and upper integrals in every
%           direction, so the parameters pair up under the index map
%           1->1, 2->3, 3->2;
%
% OUTPUT
% - prog:   Same LPI program structure as the input, but with the added
%           constraints enforcing P==0.
%
% NOTES
% To enforce P==Q, call LPI_EQ_SDOPVAR(PROG,P-Q); '@sdopvar/plus' aligns the
% monomial and decision variable bases of the two operators.
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
dvars = cellstr(string(P.Zd(:)));
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
if q>0
    known = cellstr(string(prog.decvartable(:)));
    is_new = ~ismember(dvars,known);
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

% % % Impose the constraints, one parameter at a time
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
    M = [reshape(full(A_eff),1,[]); B_eff];
    [sg,ii,vv] = find(M);
    Cdp = sparse(sg,ii,vv,q+1,n_eff);
    Dk = dpvar(Cdp,zeros(1,0),{},dvars,[1,n_eff]);
    prog = soseq(prog,Dk);
end

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
