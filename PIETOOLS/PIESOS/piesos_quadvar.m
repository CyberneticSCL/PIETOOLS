function [prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog,ZL,ZR,opts)
% [PROG,V,PMAT,ZOPL,ZOPR] = PIESOS_QUADVAR(PROG,ZL,ZR,OPTS) computes a
% distributed polynomial decision variable based on the vectors of
% distributed polynomials ZL and ZR, typically bases of distributed
% monomials, taking the quadratic form
%   V = <ZOPL*ZL, PMAT*ZOPR*ZR>
% where ZOPL and ZOPR are bases of tensor-PI operators, acting on the
% entries of ZL and ZR, and PMAT(s) is a matrix of polynomials in primary
% variables s, parameterized by decision variables.
% PMAT is a general matrix by default, a symmetric matrix if
% opts.sym==true, and an SOS matrix if opts.pos==true. For ZL=ZR and
% opts.pos==true, V is a distributed SOS polynomial.
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           decision variable;
% - ZL:     nL x 1 'polyopvar' object representing a vector of distributed
%           polynomials, e.g. a basis of distributed monomials, acted on by
%           the operators ZopL: block i of ZopL acts on entry i of ZL, see
%           NOTES (for an entry of otype 3, both output types, two block
%           rows of ZopL act on it);
% - ZR:     (optional) nR x 1 'polyopvar' object representing a vector of
%           distributed polynomials, acted on by the operators ZopR. If
%           omitted or empty, the same basis is used on the left and the
%           right, ZR=ZL;
% - opts:   (optional) struct with fields
%     pos:  set to true to declare Pmat as a positive semidefinite (SOS)
%           matrix. Requires the numbers of rows of the blocks of ZopL and
%           ZopR, and the degrees degs.L.out and degs.R.out, to match.
%           Defaults to false;
%     sym:  set to true to declare Pmat as a symmetric matrix, under the
%           same requirements. Defaults to false;
%   psatz:  (optional, only if opts.pos==true) 1 x N array of integers, or
%           logical array, specifying psatz multipliers g_k(s)>=0 on the
%           domain, with Pmat = sum_k g_k*Pmat_k and Pmat_k an SOS matrix.
%           If psatz(i)=qi, the multipliers are (si-ai)^j*(bi-si)^(qi-j),
%           j in {0,...,qi};
%    degs:  array, cell, or struct with fields 'in' and 'out', specifying
%           the degrees of the monomials in the operators ZopL and ZopR,
%           as in `piesos_quadopvar`. The degrees in degs.in define the
%           monomials in the dummy variables, in terms of which the
%           kernels of ZopL and ZopR are defined, those in degs.out the
%           monomials in the primary variables, defining the basis of
%           Pmat(s). The same degrees are used for ZopL and ZopR, unless
%           degs is a struct with fields L and R, specifying the degrees
%           of each. Defaults to 0. A field 'deg' is used as 'degs' if
%           'degs' is not specified, and ignored otherwise, with a
%           warning;
%    L/R:   (optional) structs with fields 'degs', 'exclude', 'sep', and
%           'otype', specifying separate options for the operators ZopL
%           (L) and ZopR (R): 'exclude' the multiplier, lower or upper
%           integral terms to exclude, 'sep' whether to enforce
%           separability along each variable, 'otype' the type of product
%           of the factors (1 pointwise, 2 tensor product, 3 both). The
%           degrees of the left and right basis may be specified in
%           opts.L.degs and opts.R.degs, instead of opts.degs.L and
%           opts.degs.R (specifying them for the same basis in both
%           places is an error). If only one of L and R is specified, it
%           is used for both, except that a specified opts.degs still
%           gives the degrees of the other operator. The fields
%           'exclude', 'sep', and 'otype' of opts are used for each
%           operator for which L or R does not specify them;
% exclude, sep:  (optional) options of the operators ZopL and ZopR, as in
%           `piesos_quadopvar`;
%   otype:  (optional) type of output of ZopL and ZopR: 1 for a pointwise
%           product, 2 for a tensor product, 3 for both (see `monopmat`),
%           used for both operators unless opts.L or opts.R specify it.
%           Defaults to 1;
%  intmat:  (optional) regions of integration of the copies si_1,...,si_di
%           of the output variables, a 1 x N cell of qi x 2*di arrays or
%           presets 'full' (default), 'steps', 'all', see
%           `piesos_quadopvar`. V sums a Gram form over each of the
%           q1*...*qN combinations of regions, each with its own Pmat;
%  unique:  (optional) logical, true to use the pointwise product (otype 1,
%           or the pointwise part of otype 3) for only one of the entries
%           that are monomials of the same factors in different orders,
%           e.g. z o y and y o z (`dmonomials` with 'full'), as their
%           pointwise products are the same function, see NOTES. Dropping
%           an entry of otype 1 is an error if degs, exclude or sep are
%           given per entry (as cells). Defaults to true;
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the decision
%           variables defining Pmat added (and Pmat>=0, for opts.pos);
% - V:      scalar 'polyopvar' object representing the distributed
%           polynomial functional
%               V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>;
% - Pmat:   rL x rR cell, for rL and rR the numbers of block rows of ZopL
%           and ZopR, with element Pmat{i,j} an mi x nj 'dpvar' object in
%           the primary variables, coupling the ith block row of ZopL and
%           the jth block row of ZopR, see `piesos_quadopvar`. With
%           regions other than the whole domain (opts.intmat), a
%           q1 x ... x qN cell of such cells, Pmat{k1,...,kN} for rows
%           k1,...,kN of the regions;
% - ZopL/ZopR:  'tensopmat' objects representing the bases of operators on
%           the entries of ZL and ZR (block i acts on the output of row i
%           of ZL.C);
%
% NOTES
% - Entry i of ZL is ZL_i(x) = sum_j C_ij*Z_j(x), for the monomials Z_j in
%   the rows of ZL.Z and the coefficients ZL.C of ZL (a single term for a
%   monomial basis). Block i of ZopL acts on the output of the
%   coefficients of row i, so that V is defined by the operators
%   ZopL_i*C_ij, the blocks of the product ZopL*ZL.C (mixed product), on
%   the monomials Z_j.
% - The coefficients must be tensor-PI operators with a Kronecker output
%   (otype(1)=2), so that the monomials are tensor products of the
%   variables (x1 o x1); elementwise products, as in x1.*x1, are not
%   supported. A factor that is itself a pointwise product, nested in a
%   tensor product, is not supported either. The factors of a monomial
%   must have the same size, and its states depend on the same variables.
% - The states in ZL and ZR are matched by name, and must have the same
%   size and domain in both.
% - For full generality, a basis should include distributed monomials in
%   all orders, e.g. both y o z and z o y. For pointwise products
%   (otype=1), however, y o z and z o y involve the same elements
%   (perhaps in different order): opts.unique removes these duplicates
%   for otype 1 and 3; set it to false to keep them.
% - Unrecognized fields of opts give a warning (PIESOS:unknownOption).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_quadvar
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

% Verify that the monomials are appropriately specified
is_eqZLR = false;
if nargin<=2 || isempty(ZR)
    % Assume same basis on left and right
    is_eqZLR = true;
    ZR = ZL;
end
if ~isa(ZL,'polyopvar') || ~isa(ZR,'polyopvar')
    error("Bases of distributed monomials must be specified as 'polyopvar' objects.")
end
if ~ZL.chkval
    % Guarantee Z1.Z/Z1.C are mutually consistent
    ZL = initialize(ZL);
end
if ~is_eqZLR
    % Guarantee Z2.Z/Z2.C are mutually consistent
    if ~ZR.chkval
        ZR = initialize(ZR); 
    end
    is_eqZLR = isequal(ZL,ZR);
end
% The operators ZopL, ZopR act on the output of the coefficients ZL.C, ZR.C
% (e.g. T*x): extract the function space of the output of each row, and its
% nonzero columns
CL = ZL.C;      CR = ZR.C;
[dimsL,depmatL,colsL] = coefficient_blocks(CL,'L');
if is_eqZLR
    dimsR = dimsL;      depmatR = depmatL;      colsR = colsL;
else
    [dimsR,depmatR,colsR] = coefficient_blocks(CR,'R');
end
% Declare default monomial degree if not specified
if nargin<=3
    opts = struct();
end
if ~isa(opts,'struct')
    error("Options must be specified as a 'struct'.")
end
check_opts(opts,'piesos_quadvar');
% Allow degrees to be specified using 'deg' if no 'degs' are specified
if isfield(opts,'deg')
    if isfield(opts,'degs')
        warning('PIESOS:quadvar:renamedDeg', ...
            "Option 'deg' is ignored, as opts.degs is specified.");
    else
        warning('PIESOS:quadvar:renamedDeg', ...
            "Option 'deg' is renamed 'degs': using opts.degs = opts.deg.");
        opts.degs = opts.deg;
    end
    opts = rmfield(opts,'deg');
end
% Check options for the operators
% Degrees: in opts.degs, or for the left and right basis in opts.L.degs and
% opts.R.degs, which take precedence over a general opts.degs. They are 
% passed on to piesos_quadopvar as degs, with fields L and R if specified 
% separately
has_degs = isfield(opts,'degs');
has_degsL = isfield(opts,'L') && isstruct(opts.L) && isfield(opts.L,'degs');
has_degsR = isfield(opts,'R') && isstruct(opts.R) && isfield(opts.R,'degs');
if has_degsL || has_degsR
    if has_degs && isstruct(opts.degs) && any(isfield(opts.degs,{'L','R'}))
        % Separate degrees for the left and right basis in opts.degs
        degs = opts.degs;
        if (has_degsL && isfield(degs,'L')) || (has_degsR && isfield(degs,'R'))
            error("Degrees of the left or right basis cannot be specified both in opts.degs and in opts.L or opts.R.")
        end
    elseif has_degs
        % The same degrees for the left and right basis, unless specified below
        degs = struct('L',{opts.degs},'R',{opts.degs});
    else
        degs = struct();    % if only one of the bases, the other has the same
    end
    if has_degsL
        degs.L = opts.L.degs;
    end
    if has_degsR
        degs.R = opts.R.degs;
    end
elseif has_degs
    degs = opts.degs;
else
    degs = 0;
end
% Remove the degrees from the options of the operators, and also opts.L, 
% opts.R if they only specified degrees, so that opts.exclude etc. apply to 
% both
only_degsL = has_degsL && isscalar(fieldnames(opts.L));
only_degsR = has_degsR && isscalar(fieldnames(opts.R));
if has_degs
    opts = rmfield(opts,'degs');
end
if has_degsL
    opts.L = rmfield(opts.L,'degs');
end
if has_degsR
    opts.R = rmfield(opts.R,'degs');
end
if (only_degsL || ~isfield(opts,'L')) && (only_degsR || ~isfield(opts,'R'))
    if only_degsL
        opts = rmfield(opts,'L');
    end
    if only_degsR
        opts = rmfield(opts,'R');
    end
end

% The pointwise product once for monomials of the same factors in other
% orders (opts.unique, see NOTES): an entry whose factors appear in an
% earlier entry gets otype 2 instead of 3, or is left out for otype 1
if ~isfield(opts,'unique') || opts.unique
    % A lone opts.L or opts.R applies to both bases (as in piesos_quadopvar),
    % also for the otype of each side here
    if isfield(opts,'L') && ~isfield(opts,'R')
        opts.R = opts.L;
    elseif isfield(opts,'R') && ~isfield(opts,'L')
        opts.L = opts.R;
    end
    for side = 'LR'
        if side=='L'
            Zs = ZL;    cols = colsL;
        else
            Zs = ZR;    cols = colsR;
        end
        otype = 1;
        if isfield(opts,'otype')
            otype = opts.otype;
        end
        if isfield(opts,side) && isfield(opts.(side),'otype')
            otype = opts.(side).otype;
        end
        otype = otype(:,1);
        if isscalar(otype)
            otype = repmat(otype,numel(cols),1);
        end
        [otype_new,is_kept] = pointwise_once(Zs.Z,Zs.C,cols,otype);
        if isequal(otype_new,otype) && all(is_kept)
            continue        % no entries of the same factors in other orders
        end
        if ~all(is_kept)
            % Options given per entry would no longer match the entries: a
            % cell with an element per entry, as degs, degs.in/out, or within
            % degs.L/R (a 1 x 1 cell is used for all entries)
            per_entry = @(D) iscell(D) && numel(D)>1;
            has_cells = per_entry(degs);
            if isstruct(degs)
                for fn = fieldnames(degs)'
                    Df = degs.(fn{1});
                    has_cells = has_cells || per_entry(Df) || ...
                        (isstruct(Df) && any(cellfun(@(g) per_entry(Df.(g)),fieldnames(Df))));
                end
            end
            for fname = {'exclude','sep'}
                has_cells = has_cells || (isfield(opts,fname{1}) && iscell(opts.(fname{1}))) || ...
                    (isfield(opts,side) && isfield(opts.(side),fname{1}) && iscell(opts.(side).(fname{1})));
            end
            if has_cells
                error("An entry of otype 1 with the same factors as another (opts.unique) cannot be left out with options given per entry: set opts.unique=false, or remove it from the basis.")
            end
            % Block rows ({}): '()' would index inside the coefficients, which
            % fails for Kronecker products, e.g. with an ODE factor
            Zs = Zs{find(is_kept),:};   % numeric indices: polyopvar      %#ok<FNDSB>
            [dims_s,depmat_s,cols] = coefficient_blocks(Zs.C,side);
            if side=='L'
                ZL = Zs;    CL = Zs.C;  dimsL = dims_s;     depmatL = depmat_s;     colsL = cols;
            else
                ZR = Zs;    CR = Zs.C;  dimsR = dims_s;     depmatR = depmat_s;     colsR = cols;
            end
        end
        % Both opts.L and opts.R, as piesos_quadopvar uses one given alone
        % for both, and fills fields not given from opts
        if isfield(opts,'L') && ~isfield(opts,'R')
            opts.R = opts.L;
        elseif isfield(opts,'R') && ~isfield(opts,'L')
            opts.L = opts.R;
        elseif ~isfield(opts,'L')
            opts.L = struct();
            opts.R = struct();
        end
        opts.(side).otype = otype_new(is_kept);
    end
    is_eqZLR = is_eqZLR && isequal(ZL,ZR);
end

% Extract the domains. The operators act on the output of the coefficients,
% the output spaces of the rows of CL, CR
dims = struct('R',dimsR,'L',dimsL);
vars = struct('R',{CR.vars.out},'L',{CL.vars.out});
dom = struct('R',CR.dom.out,'L',CL.dom.out);
depmat = struct('R',depmatR,'L',depmatL);
if any(arrayfun(@(i) numel(unique(nonzeros(depmat.L(i,:))))>1,1:size(depmat.L,1))) || ...
        any(arrayfun(@(i) numel(unique(nonzeros(depmat.R(i,:))))>1,1:size(depmat.R,1)))
    error("Distributed monomials of states depending on different variables are not supported.")
end

% Declare a basis of operators on this domain, and a matrix Pmat(s) in the
% primary variables (nL x nR cell, a block row for each entry of ZL and
% a block column for each of ZR), so that Pop = ZopL'*Pmat*ZopR; or one
% Pmat for each region of opts.intmat, Pregs{k}
[prog,Pmat,ZopL,ZopR,intmat] = piesos_quadopvar(prog,dims,vars,dom,depmat,degs,opts);
Pregs = Pmat;
if ~iscell(Pmat{1})
    Pregs = {Pmat};     % the whole domain
end

% The terms of the pairs of entries (ridx,cidx) and (cidx,ridx) are equal
% if Pmat is symmetric and the same entries and operators appear on the
% left and right, in which case only the upper triangle is computed
is_sym = (isfield(opts,'pos') && opts.pos) || (isfield(opts,'sym') && opts.sym);
is_sym = is_sym && is_eqZLR && isequal(ZopL,ZopR);
% The block rows of ZopL, ZopR (of Pmat), and the entry of ZL, ZR on
% which each acts (its only nonempty block): entry i for block row i, or
% two block rows for an entry of otype 3 (see piesos_quadopvar)
[nZL,nZR] = size(Pregs{1});
[~,blk_idcsL] = max(~cellfun(@isempty,ZopL.ops),[],2);
[~,blk_idcsR] = max(~cellfun(@isempty,ZopR.ops),[],2);
if is_sym
    [ridcs,cidcs] = find(triu(true(nZL,nZR)));
else
    [ridcs,cidcs] = find(true(nZL,nZR));
end
% Express the monomials in ZL and ZR in terms of the same states
if is_eqZLR
    ZLc = ZL;       ZRc = ZR;
else
    [ZLc,ZRc] = common_states(ZL,ZR);
end

% Compose the basis operators with the coefficients by one product: block (i,j)
% of ZopL*ZL.C is ZopL_i*C_ij, acting on the monomial Z_j by the mixed product of
% the factors. The empty parameters of the excluded terms are filled first, as
% the 'sopvar' product fails on them
PL = fill_empty_params(ZopL)*CL;
if is_eqZLR && isequal(ZopL,ZopR)
    PR = PL;
else
    PR = fill_empty_params(ZopR)*CR;
end

% A term for each pair of rows (ridx,cidx) of Pmat and each pair of nonzero blocks
% (a,b) in these rows, acting on the monomial [ZL.Z(a,:),ZR.Z(b,:)]: the columns
% of terms are ridx, cidx, a, b
terms = zeros(0,4);
for k=1:numel(ridcs)
    % A block row without rows (all its terms excluded by opts.exclude, e.g. the
    % multiplier of a state without variables) contributes no term
    if size(ZopL.ops{ridcs(k),blk_idcsL(ridcs(k))},1)==0 || size(ZopR.ops{cidcs(k),blk_idcsR(cidcs(k))},1)==0
        continue
    end
    [a,b] = ndgrid(colsL{blk_idcsL(ridcs(k))},colsR{blk_idcsR(cidcs(k))});
    terms = [terms; repmat([ridcs(k),cidcs(k)],numel(a),1),a(:),b(:)];
end
if isempty(terms)
    error("All terms of the functional are excluded by the option 'exclude'.")
end

% Compute the list of monomials that appear in the distributed polynomial.
Z_full = [ZLc.Z(terms(:,3),:), ZRc.Z(terms(:,4),:)];
for i=1:size(Z_full,1)
    Z_i = Z_full(i,:);
    Z_full(i,:) = [Z_i(Z_i>0), Z_i(Z_i==0)];
end
[M,Z_lin] = uniquerows_integerTable([sum(Z_full>0,2),Z_full]);  % Z_full = M*Z_lin
Z_lin = Z_lin(:,2:end);
old2new_idcs = M*(1:size(Z_lin,1))';               % Z_full = Z_lin(old2new_idcs,:);
nZ_lin = size(Z_lin,1);

% Initialize an empty coefficient operator acting on the new basis
Kcell = num2cell(zeros(1,nZ_lin));

% Regions: for Pregs{k}, a row of intmat{i} for each output variable i,
% the region of all its copies over which V integrates
nreg = [cellfun(@(A) size(A,1),intmat),1];
Q = numel(Pregs);
intmat_reg = cell(1,Q);     rows = cell(1,numel(nreg));
for k=1:Q
    [rows{:}] = ind2sub([nreg,1],k);  % nreg=1 without output variables
    intmat_reg{k} = cellfun(@(A,r) A(r,:),intmat,rows(1:end-1),'UniformOutput',false);
end

% Determine the operator acting on each of the distributed monomials
% One iteration for each term (ridx,cidx,a,b), with the weight Wij of a row pair
% computed once for its terms
for t=1:size(terms,1)
    ridx = terms(t,1);      cidx = terms(t,2);
    if t==1 || ~isequal(terms(t,1:2),terms(t-1,1:2))
        % Convert the block of decision variables to the linear format: Kronecker
        % (pointwise) product along the output dim, proper tensor product along
        % the input dim. The factors are ordered as the monomial [ZL(a),ZR(b)] in
        % Z_full. For symmetric Pmat, (cidx,ridx) and (ridx,cidx) give the same term
        % m, n: row counts of the block rows ridx of ZopL and cidx of ZopR
        m = size(ZopL.ops{ridx,blk_idcsL(ridx)},1);
        n = size(ZopR.ops{cidx,blk_idcsR(cidx)},1);
        % Pij(1:end) as a row, for each region
        Wij = cell(1,Q);
        for k=1:Q
            Wij{k} = Pregs{k}{ridx,cidx}(:).';
            if is_sym && cidx>ridx
                Wij{k} = 2*Wij{k};
            end
        end
    end
    Cop = otimes(repmat(PL.ops{ridx,terms(t,3)},n,1),repelem(PR.ops{cidx,terms(t,4)},m,1),[1,2]);

    % Add the terms of all regions to the correct distributed monomial
    mon_idx = old2new_idcs(t);
    for k=1:Q
        Kcell{mon_idx} = Kcell{mon_idx} + int(Cop,Wij{k},intmat_reg{k});
    end
end

% Declare the coefficient matrix defining V
Kop = tensopmat(Kcell);

% Declare the polynomial defined by coefficients K and monomials Z
V = ZLc;
V.Z = Z_lin;
% The output space of V is set by the assignment of V.C
V.C = Kop;
V.chkval = true;

% Get rid of duplicate terms in the functional
V = combine_terms(V);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [dims,depmat,cols] = coefficient_blocks(C,side)
% [DIMS,DEPMAT,COLS] = COEFFICIENT_BLOCKS(C,SIDE) checks the coefficients C
% of a vector of distributed polynomials (tensor-PI operators with a
% Kronecker output) and returns the function space of the output of each
% row of C, on which the basis operator of the row acts, and its nonzero
% columns: the size DIMS(i) of each factor of the output of row i (the size
% of the output, if its factors cannot be distinguished, or it has no
% spatial variables), the number DEPMAT(i,:) of copies of each output
% variable, and the indices COLS{i} of the monomials in row i. SIDE ('L' or
% 'R') is used in error messages.
%
% INPUTS
% - C:      nr x nZ 'tensopmat', the coefficients of the vector;
% - side:   char, 'L' or 'R';
%
% OUTPUTS
% - dims:   nr x 1 array, the size of each factor of the output of row i;
% - depmat: nr x M array, the number of copies of each output variable in
%           the output of row i;
% - cols:   nr x 1 cell; cols{i} lists the nonzero columns of row i.

nr = size(C.ops,1);
dims = zeros(nr,1);     depmat = double(C.depmat.out);      cols = cell(nr,1);
for i=1:nr
    cols{i} = find(~cellfun(@isempty,C.ops(i,:)));
    if isempty(cols{i})
        error("Row %d of Z%s.C is zero: the entries of Z%s must be nonzero.",i,side,side)
    end
    for j=cols{i}
        blk = C.ops{i,j};
        if ~isa(blk,'tensopvar')
            error("The coefficients of Z%s must be tensor-PI operators (row %d, column %d is not).",side,i,j)
        elseif any(cellfun(@(f) isa(f,'tensopvar'),blk.ops(:)))
            error("Coefficients of Z%s with nested factors, e.g. a pointwise product of states in a tensor product, are not supported (row %d, column %d).",side,i,j)
        elseif blk.otype(1)==1 && size(blk.ops,2)>1
            error("The monomial basis Z%s must be defined in terms of tensor products of the variables, as in dmonomials, not elementwise products (row %d, column %d).",side,i,j)
        end
    end
    % The size of each factor, if each of the factors has a copy of a variable
    blk = C.ops{i,cols{i}(1)};
    dims(i) = size(blk,1);
    if any(depmat(i,:)) && size(blk.ops,2)==max(depmat(i,:))
        fdims = blk.dims(:,1);
        if any(fdims~=fdims(1))
            error("The factors of a monomial in Z%s must have the same size (row %d).",side,i)
        end
        dims(i) = fdims(1);
    elseif any(depmat(i,:)) && any(~any(blk.depmat.in,2) & ~any(blk.depmat.out,2))
        % Some factors are states without spatial variables (ODE states), 
        % which are merged into a neighboring factor with variables.
        % The size of that factor is the product of the sizes
        is_R = ~any(blk.depmat.in,2) & ~any(blk.depmat.out,2);
        fdims = blk.dims(:,1);
        fdims_eff = zeros(0,1);     acc = 1;
        for k=1:numel(fdims)
            acc = acc*fdims(k);
            if ~is_R(k)
                fdims_eff(end+1,1) = acc;       acc = 1;
            end
        end
        fdims_eff(end) = fdims_eff(end)*acc;
        if any(fdims_eff~=fdims_eff(1))
            error("The factors of a monomial in Z%s must have the same size, after merging the states without spatial variables into their neighbors (row %d).",side,i)
        end
        dims(i) = fdims_eff(1);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function T = fill_empty_params(T)
% T = FILL_EMPTY_PARAMS(T) replaces the empty parameters params{k}=[] of
% the 'sopvar' objects in the tensor-PI operator or the 'tensopmat' T
% (block diagonal, or with the blocks of both output types stacked), left
% empty for the excluded multiplier, lower and upper terms, by zero
% matrices of the right size.
%
% INPUTS
% - T:  'tensopmat', 'tensopvar' or 'sopvar' object;
%
% OUTPUTS
% - T:  the same operator, with every empty parameter of its 'sopvar'
%       factors a sparse zero matrix.

% Containers: fill the parameters of each nonempty element
if isa(T,'tensopmat')
    for k=find(~cellfun(@isempty,T.ops(:)))'
        T.ops{k} = fill_empty_params(T.ops{k});
    end
    T.chkval = true;
elseif isa(T,'tensopvar')
    for k=1:numel(T.ops)
        T.ops{k} = fill_empty_params(T.ops{k});
    end
% 'sopvar': empty parameters become zero matrices of
% dims(1)*(number of left monomials) x dims(2)*(number of right monomials)
elseif isa(T,'sopvar')
    params = T.params;
    nL = prod(cellfun(@numel,T.ZL));    nR = prod(cellfun(@numel,T.ZR));
    for k=1:numel(params)
        if isempty(params{k})
            params{k} = sparse(T.dims(1)*nL,T.dims(2)*nR);
        end
    end
    T = sopvar(params,T.vars,T.ZL,T.ZR,T.dom,T.dims);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [otype,is_kept] = pointwise_once(Zmat,C,cols,otype)
% [OTYPE,IS_KEPT] = POINTWISE_ONCE(ZMAT,C,COLS,OTYPE) keeps the pointwise
% product (otype 1, or the pointwise part of otype 3) of only one of the
% entries of a basis that are monomials of the same factors in different
% orders, e.g. z o y and y o z, whose pointwise products are the same
% function (opts.unique of piesos_quadvar). The entry in the order of the
% states keeps it, else the first.
%
% INPUTS
% - Zmat:   nZ x dmax array, row j the states of the factors of monomial j
%           (zeros after the last factor), the field Z of the basis;
% - C:      nr x nZ 'tensopmat', the coefficients of the basis;
% - cols:   nr x 1 cell, cols{i} the monomials of entry i (see
%           coefficient_blocks);
% - otype:  nr x 1 array, the type of output of each entry: 1 pointwise,
%           2 tensor product, 3 both;
%
% OUTPUTS
% - otype:  nr x 1 array, with 3 replaced by 2 for an entry whose pointwise
%           product is that of an entry kept before it;
% - is_kept: nr x 1 logical, false for such an entry of otype 1, which has
%           no rows left;

nr = numel(cols);
% The states of an entry that is a single monomial with identity
% coefficient and two or more factors, sorted; empty for other entries
fctrs_sorted = cell(nr,1);
for i=1:nr
    if isscalar(cols{i})
        fctrs = nonzeros(Zmat(cols{i},:)).';
        if numel(fctrs)>1 && is_id_op(C.ops{i,cols{i}})
            fctrs_sorted{i} = sort(fctrs);
        end
    end
end
is_kept = true(nr,1);
is_done = cellfun(@isempty,fctrs_sorted);
for i=find(~is_done).'
    if is_done(i)
        continue
    end
    % The entries of the same factors, those in the order of the states first
    grp = find(cellfun(@(f) isequal(f,fctrs_sorted{i}),fctrs_sorted));
    is_ordered = arrayfun(@(k) issorted(nonzeros(Zmat(cols{k},:))),grp);
    grp = [grp(is_ordered);grp(~is_ordered)];
    has_pointwise = false;
    for k=grp.'
        if has_pointwise && otype(k)==3
            otype(k) = 2;
        elseif has_pointwise && otype(k)==1
            is_kept(k) = false;
        elseif otype(k)~=2
            has_pointwise = true;
        end
    end
    is_done(grp) = true;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %