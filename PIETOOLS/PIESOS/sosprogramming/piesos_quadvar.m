function [prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog,ZL,ZR,opts)
% [PROG,V,PMAT,ZOPL,ZOPR] = PIESOS_QUADVAR(PROG,ZL,ZR,OPTS) computes a
% distributed polynomial decision variable based on the vectors of
% distributed polynomials ZL and ZR, typically bases of distributed
% monomials, taking the quadratic form
%   V = <ZOPL*ZL, PMAT*ZOPR*ZR>
% where ZOPL and ZOPR are bases of tensor-PI operators, acting on the
% entries of ZL and ZR, and PMAT(s) is a matrix of
% polynomials in primary variables s, parameterized by decision variables.
% PMAT is a general matrix by default, a symmetric matrix if
% opts.sym==true, and an SOS matrix if opts.pos==true. For ZL=ZR and
% opts.pos==true, V is a distributed SOS polynomial.
%
% INPUTS
% - prog:   PIESOS program structure in which to declare the decision
%           variable;
% - ZL:     nL x 1 'polyopvar' object representing a vector of distributed
%           polynomials, e.g. a basis of distributed monomials, acted on by
%           the operators ZopL: block i of ZopL acts on entry i of ZL, see
%           NOTES;
% - ZR:     (optional) nR x 1 'polyopvar' object representing a vector of
%           distributed polynomials, acted on by the operators ZopR. If
%           omitted or empty, the same basis is used on the left and the
%           right, ZR=ZL;
% - opts:   (optional) struct with fields
%     pos:  set to true to declare Pmat as a positive semidefinite (SOS)
%           matrix. Requires the numbers of rows of the blocks of ZopL and
%           ZopR, and the degrees degs.L.out and degs.R.out, to match, see
%           `piesos_quadopvar`. Defaults to false;
%     sym:  set to true to declare Pmat as a symmetric matrix, under the
%           same requirements. Defaults to false;
%   psatz:  (optional, only if opts.pos==true) 1 x N array of integers, or
%           logical array, specifying psatz multipliers g_k(s)>=0 on the
%           domain, with Pmat = sum_k g_k*Pmat_k and Pmat_k an SOS matrix.
%           If psatz(i)=qi, the multipliers are (si-ai)^j*(bi-si)^(qi-j),
%           j in {0,...,qi}, see `piesos_quadopvar`;
%    degs:  array, cell, or struct with fields 'in' and 'out', specifying
%           the degrees of the monomials in the operators ZopL and ZopR,
%           as in `piesos_quadopvar`. The degrees in degs.in define the
%           monomials in the dummy variables, in terms of which the kernels
%           of ZopL and ZopR are defined, those in degs.out the monomials 
%           in the primary variables, defining the basis of Pmat(s). The 
%           same degrees are used for ZopL and ZopR, unless degs is a 
%           struct with fields L and R, specifying the degrees of each;
%    L/R:   (optional) structs with fields 'degs', 'exclude', 'sep', and
%           'otype', specifying separate options for the operators ZopL
%           (L) and ZopR (R), see `monopmat` for the latter three. The
%           degrees of the left and right basis may be specified in
%           opts.L.degs and opts.R.degs, instead of opts.degs.L and
%           opts.degs.R (specifying them for the same basis in both
%           places is an error). If only one of L and R is specified, it
%           is used for both. If neither is specified, or only to
%           specify degrees, the fields 'exclude', 'sep', and 'otype' of
%           opts are used for both operators;
%
% OUTPUTS
% - prog:   PIESOS program structure with the decision variables defining
%           Pmat added;
% - V:      scalar 'polyopvar' object representing the distributed
%           polynomial functional
%               V(x) = <ZopL*ZL(x), Pmat*ZopR*ZR(x)>;
% - Pmat:   nL x nR cell, with element Pmat{i,j} an mi x nj 'dpvar' object
%           in the primary variables, coupling the ith block row of ZopL
%           and the jth block row of ZopR, i.e. entry i of ZL and entry j
%           of ZR, see `piesos_quadopvar`;
% - ZopL/ZopR:  'tensopmat' objects representing the bases of operators on
%           the entries of ZL and ZR (block i acts on the output of row i
%           of ZL.C), see `piesos_quadopvar`;
%
% NOTES
% - Entry i of ZL is ZL_i(x) = sum_j C_ij*Z_j(x), for the monomials Z_j in
%   the rows of ZL.Z and the coefficients ZL.C of ZL (a single term for a
%   monomial basis). Block i of ZopL acts on the output of the coefficients
%   of row i, so that V is defined by the operators ZopL_i*C_ij, the blocks
%   of the product ZopL*ZL.C (mixed product), on the monomials Z_j. 
% - The coefficients must be tensor-PI operators with a Kronecker output
%   (otype(1)=2), so that the monomials are tensor products of the
%   variables, as in dmonomials; elementwise products, as in x1.*x1, are
%   not supported. A factor that is itself a pointwise product, nested in a
%   tensor product, is not supported either. The factors of a monomial must
%   have the same size, and its states depend on the same variables.
% - The states in ZL and ZR are matched by name, and must have the same
%   size and domain in both.
% - If Pmat is symmetric, ZL=ZR, and ZopL=ZopR, only the upper triangle of 
%   Pmat is used to compute V. 

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
% DJ, 10/15/2026: Initial coding;
% DJ, 10/22/2026: Support coefficients ZL.C, ZR.C, e.g. T*x;

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
% a block column for each of ZR), so that Pop = ZopL'*Pmat*ZopR;
[prog,Pmat,ZopL,ZopR] = piesos_quadopvar(prog,dims,vars,dom,depmat,degs,opts);

% The terms of the pairs of entries (ridx,cidx) and (cidx,ridx) are equal
% if Pmat is symmetric and the same entries and operators appear on the
% left and right, in which case only the upper triangle is computed
is_sym = (isfield(opts,'pos') && opts.pos) || (isfield(opts,'sym') && opts.sym);
is_sym = is_sym && is_eqZLR && isequal(ZopL,ZopR);
nZL = numel(dimsL);     nZR = numel(dimsR);
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
% @sopvar/mtimes fails on them
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
    [a,b] = ndgrid(colsL{ridcs(k)},colsR{cidcs(k)});
    terms = [terms; repmat([ridcs(k),cidcs(k)],numel(a),1),a(:),b(:)];
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
        m = size(ZopL.ops{ridx,ridx},1);    n = size(ZopR.ops{cidx,cidx},1);    % row counts of ZopL_i, ZopR_j
        Wij = vec_dpvar(Pmat{ridx,cidx});       % Pij(1:end) as a row
        if is_sym && cidx>ridx
            Wij = 2*Wij;
        end
    end
    Kij = int(otimes(repmat(PL.ops{ridx,terms(t,3)},n,1),repelem(PR.ops{cidx,terms(t,4)},m,1),[1,2]),Wij);

    % Match the obtained kernel with the correct distributed monomial
    mon_idx = old2new_idcs(t);
    Kcell{mon_idx} = Kcell{mon_idx} + Kij;
end

% Declare the coefficient matrix defining V
Kop = tensopmat(Kcell);

% Declare the polynomial defined by coefficients K and monomials Z
V = ZLc;
V.Z = Z_lin;
V.C = Kop;
V.dims.out = 1;
V.vars.out = Kop.vars.out;
V.dom.out = Kop.dom.out;
V.depmat.out = Kop.depmat.out;
V.chkval = true;

% Get rid of duplicate terms in the functional
V = combine_terms(V);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% 
function W = vec_dpvar(P)
% W = VEC_DPVAR(P) returns the 1 x numel(P) 'dpvar' W such that W = P(:).'
% in the order of the linear index (column-major), computed from the
% coefficients of P. For P of dimension m x n with q decision variables, 
% nZ monomials, P.C is m*(q+1) x n*nZ, and W.C is (q+1) x m*n*nZ, 
% entry (i,j) of P moving to block k = i+(j-1)*m of columns.
if ~isa(P,'dpvar')
    W = P(1:end);
    if size(W,1)>1,     W = W.';    end
    return
end
[m,n] = size(P);
nd = numel(P.dvarname)+1;       nZ = size(P.degmat,1);
[r,c,v] = find(P.C);
i = ceil(r/nd);                 d = r-(i-1)*nd;
j = ceil(c/nZ);                 z = c-(j-1)*nZ;
k = i+(j-1)*m;
W = compress(dpvar(sparse(d,(k-1)*nZ+z,v,nd,m*n*nZ),P.degmat,P.varname,P.dvarname,[1,m*n]));   % as dpvar/subsref
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
        elseif size(blk.ops,2)>1 && ~any(depmat(i,:))
            % The factors of the block cannot be composed with the single factor of ZopL_i
            error("Products of states without spatial variables in the basis Z%s are not supported (row %d, column %d).",side,i,j)
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
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function T = fill_empty_params(T)
% T = FILL_EMPTY_PARAMS(T) replaces the empty parameters params{k}=[] of
% the 'sopvar' objects in the tensor-PI operator or the diagonal 
% 'tensopmat' T, left by monopmat for the excluded multiplier, lower and 
% upper terms, by zero matrices of the right size.

if isa(T,'tensopmat')
    for i=1:size(T.ops,1)
        T.ops{i,i} = fill_empty_params(T.ops{i,i});
    end
    T.chkval = true;
elseif isa(T,'tensopvar')
    for k=1:numel(T.ops)
        T.ops{k} = fill_empty_params(T.ops{k});
    end
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