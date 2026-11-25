function [prog,Pmat,ZopL,ZopR,intmat] = piesos_quadopvar(prog,dims,vars,dom,depmat,degs,opts)
% [PROG,PMAT,ZOPL,ZOPR,INTMAT] = PIESOS_QUADOPVAR(PROG,DIMS,VARS,DOM,
% DEPMAT,DEGS,OPTS) declares an operator variable of the form
%   Pop = ZOPL'*PMAT*ZOPR
% where PMAT(s) is a polynomial or SOS function in primary variables s,
% and ZOPL and ZOPR are bases of tensor-PI operators, acting on the
% function spaces
%   L2^{m1}[OmegaL_1] x ... x L2^{mp}[OmegaL_p]
%   L2^{n1}[OmegaR_1] x ... x L2^{nq}[OmegaR_q]
% defined by DIMS, VARS, DOM, and DEPMAT. (The operator Pop itself is not
% returned, since it cannot be represented by a 'tensopmat'.)
%
% INPUTS
% - prog:   'struct' specifying a PIESOS program in which to declare the
%           operator variable;
% - dims:   struct with fields
%  L/out:   p x 1 array specifying the number of coupled state variables
%           (mi) on which each of the blocks in ZopL acts;
%   R/in:   q x 1 array specifying the number of coupled state variables
%           (ni) on which each of the blocks in ZopR acts;
% - vars:   struct with fields
%  L/out:   1 x M cellstr specifying the names of the spatial variables
%           in terms of which the outputs are defined;
%   R/in:   1 x N cellstr specifying the names of the spatial variables
%           in terms of which the inputs are defined;
% - dom:    struct with fields
%  L/out:   M x 2 array specifying the intervals on which the output
%           spatial variables are defined;
%   R/in:   N x 2 array specifying the intervals on which the input spatial
%           variables are defined;
% - depmat: struct with fields
%  L/out:   p x M integer array specifying the degree to which each of the
%           states depends on each of the output variables, i.e. the
%           distributed monomial degree of the state;
%   R/in:   q x N integer array specifying the degree to which each of the
%           states depends on each of the input variables (d_in), i.e. the
%           distributed monomial degree of the state;
% - degs:   (optional) struct with fields L and R, specifying the monomial
%           degrees defining the operators ZopL and ZopR. Each field is
%           itself a struct with the following fields:
%   L/R.in: 1 x (p or q) cell, specifying the degrees of the input
%           variables in each operator block, as degs.in of `monopmat`.
%           Element in{i} can be a di x (M or N) cell, for di the number
%           of factors, max(depmat.L/R(i,:)), specifying the monomial
%           degrees (including 0) in each factor of each input variable
%           in block i, a single row being used for all factors. It can
%           also be a di x (M or N) array of maximal monomial degrees,
%           where a single row (or a scalar, for all variables) gives the
%           maximal degrees summed over the factors, split over them;
%  L/R.out: 1 x (p or q) cell, specifying the degrees of the primary
%           variables in each operator block. Each element degs.out{i}
%           can be a 1 x P cell of column vectors, specifying the
%           exponents of the monomial basis in each of the P primary
%           variables, or a 1 x P array, specifying the maximal monomial
%           degree in each variable (a scalar is used for all variables).
%           The degrees in degs.in define the monomials in the dummy
%           variables, in terms of which the kernels of ZopL/ZopR are
%           defined. The degrees in degs.out define the monomials in the
%           primary variables, defining the basis of the SOS matrix Pmat,
%           see NOTES. If degs is an array or cell, or a struct with
%           fields 'in' and 'out' only, the same degrees are used for
%           ZopL and ZopR. For otype(1)=2, the primary variables are the
%           copies of the output variables, see NOTES. Defaults to 0;
% - opts:   (optional) struct with fields
%     pos:  set to true to declare Pmat as a positive semidefinite (SOS)
%           variable. Requires the output dimensions of ZopL and ZopR
%           (which depend on the monomial degrees), and the degrees
%           degs.L.out and degs.R.out of the primary variables, to match.
%           Defaults to false;
%     sym:  set to true to declare Pmat as a symmetric variable
%           (Pmat(s)'=Pmat(s)), under the same requirements. Defaults to
%           false;
%   psatz:  (only if opts.pos==true, ignored otherwise) 1 x P array of
%           integers, specifying the maximal order of the psatz
%           multiplier in each of the P primary variables (0 in a
%           variable on which no state depends). Here, psatz(i)=qi
%           implies inclusion of all multipliers
%               gj = (bi-si)^{qi-j}*(si-ai)^{j},    j in {0,...,qi}
%           This means a total number of prod(psatz+1) terms. Can also
%           be set to an n1 x ... x nP logical array (not for
%           otype(1)=2), including only the multipliers associated with
%           true elements: element k along dimension i corresponds to gj
%           with j=k-1 and qi=ni-1, so psatz(i)=qi is equivalent to a
%           logical array of size qi+1 along dimension i. Terms are
%           ordered by linear index in psatz. A scalar psatz is used for
%           all variables on which some state depends. For gj of order
%           qi, the degree in si of the basis monomials of Pmat is
%           reduced by ceil(qi/2), so that the degree of Pmat in si does
%           not exceed twice its degree in degs.out. On a region other
%           than the whole domain (opts.intmat), the multipliers are
%           instead, for each group of k copies related by the limits of
%           the region, the products of total order k*qi of the faces of
%           the group, e.g. of s_2-a, s_1-s_2 and b-s_1 on the triangle
%           s_2<=s_1;
% exclude, sep, otype:  options of the operators ZopL and ZopR, passed to
%           `monopmat`, used for each operator for which the option is
%           not specified in opts.L or opts.R: 'exclude' a 3^N3 logical
%           array (or a 3 x 1 array, or a 1 x q cell, one per block) of
%           the multiplier, lower and upper integral terms to exclude
%           along the N3 variables of both input and output; 'sep' a
%           1 x N3 logical array (or a 1 x q cell), true to enforce
%           separability along a variable; 'otype' the type of product of
%           the factors, a scalar or one per block, 1 by default. The
%           value otype=2 (in otype(1)) declares a tensor-product output,
%           and otype=3 both a pointwise and a tensor-product output, see
%           NOTES;
%    L/R:   (optional) structs with fields 'exclude', 'sep', and 'otype',
%           specifying separate options for the operators ZopL and ZopR.
%           If only one of L, R is specified, it is used for both. Fields
%           that are not specified are taken from opts;
%  intmat:  (optional) regions of integration of Pmat: a 1 x N cell, for
%           the N primary variables si (see NOTES; 'full' for a variable
%           on which no state depends), with element i a qi x 2*di array,
%           for di the largest number of copies of si in a block row, a
%           row per region of the copies si_1,...,si_di, as the field
%           intmat of an 'intop' (0 for ai, di+1 for bi, j for si_j, [j,j]
%           for a Dirac delta at si_j), the limits of si_j only ai, bi and
%           si_1,...,si_(j-1). A Pmat is declared for
%           each of the q1*...*qN combinations of rows. Element i may
%           also be a preset: 'full' (the whole domain), 'steps' (each
%           strict ordering of the copies), or 'all' (each ordering, with
%           ties, except all copies tied, the pointwise product of
%           otype 1). A single preset or array is used for all variables.
%           Defaults to 'full';
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the matrix
%           variable Pmat added;
% - Pmat:   rL x rR cell, for rL and rR the numbers of block rows of ZopL
%           and ZopR (p and q, or more with otype 3; rL=rR if opts.sym or
%           opts.pos), with element Pmat{i,j} an mi x nj 'dpvar' object
%           in the primary variables, where mi and nj are the numbers of
%           rows of the ith block row of ZopL.ops and the jth block row
%           of ZopR.ops. The cell represents the matrix Pmat(s)
%           parameterizing Pop = ZopL'*Pmat*ZopR, with Pmat{i,j} of the
%           form ZsL'*Q*ZsR for bases ZsL, ZsR of monomials in the
%           primary variables, and Q general, symmetric (opts.sym) or
%           positive semidefinite (opts.pos, an SOS matrix). With psatz
%           multipliers g_k, Pmat = sum_k g_k*Pmat_k, where Pmat_k is of
%           the form Zs'*Qk*Zs, with Qk >= 0 and Zs a basis of monomials
%           in the primary variables. With regions other than the whole
%           domain (opts.intmat), a q1 x ... x qN cell of such cells, with
%           Pmat{k1,...,kN} for rows k1,...,kN of intmat, each with its
%           own decision variables (and psatz terms);
% - ZopL/ZopR:  rtot x (m or n) 'tensopmat' object representing a basis of
%           operators of specified degrees on the declared function space,
%           where m/n = sum_i dims.L/R(i)^di, for di the number of factors
%           of block i (sum(dims.L/R) if every di=1), and rtot depends on
%           the monomial degrees. The kernels of the operators depend
%           only on the dummy variables;
% - intmat: 1 x N cell, the regions of opts.intmat, with presets replaced
%           by their rows, for the output variables of Zop
%           (ZopL.vars.out: those on which some state depends, sorted by
%           name), the order of the dimensions of a cell Pmat;
%
% NOTES
% - For dims, vars, dom, depmat and degs, if specified as an array, or if
%   only one of the fields L, R is specified, the same value is used for
%   ZopL and ZopR.
% - For states depending on different variables, the operators ZopL and
%   ZopR map to a function space in all of the variables appearing in
%   any state of ZopL or ZopR (the primary variables). For example, for
%   state x(s1) in ZopL and y(s1,s2) in ZopR, both map to functions of s1
%   and s2. The primary variables are ordered as vars.L, followed by those
%   of vars.R not in vars.L; degs.L/R.out, psatz and opts.intmat refer to
%   this order.
% - With otype(1)=2, Zop has a tensor-product (Kronecker) output: a block
%   acting on a monomial of d factors maps to d copies of each primary
%   variable si, named si_1,...,si_d, and V=<Zop*Z,Pmat*Zop*Z>
%   integrates over all copies. The field degs.out gives the degrees of
%   the basis Zs in each copy, as a scalar, a 1 x P array (for all
%   copies of a variable), or a d x P array (a row per copy, for the
%   largest d: blocks with fewer factors use the first rows), and a
%   psatz multiplier is used in each copy, with the same order for each
%   copy of a variable (numeric psatz only). Specify degs.in per factor,
%   e.g. [d;d] for N=1: a scalar is a total degree, split over the
%   factors.
% - Variables of different otype may be used to represent different types
%   of products. For example, int x^4 ds can only be represented using an
%   operator variable of otype=1, whereas (int x^2 ds)^2 can only be
%   represented using one of otype=2. For full generality, therefore,
%   both types of operators may need to be included in parameterizing
%   distributed (SOS) functionals.
% - With otype=3 for a block, Zop has a pointwise and a tensor-product
%   block row acting on that block (one block row if the block has a
%   single factor, for which the two coincide), so that Pmat also has
%   cross terms between the two types of products. ZopL and ZopR may
%   also have a different otype.
% - The presets of opts.intmat contain no region twice. For a basis that
%   is symmetric in the copies (of a single state, with the same degrees
%   in each factor), regions that differ only by a permutation of the
%   copies other than si_1 give the same terms, so fewer may be specified.
% - Unrecognized fields of opts give a warning (PIESOS:unknownOption).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_quadopvar
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

% Check that the inputs are appropriate
% ZopL acts on the space (L.dims,L.vars,L.dom,L.depmat), ZopR on that of R
L = struct();       R = struct();
[L.dims,R.dims] = split_LR(dims,'dims');
[L.vars,R.vars] = split_LR(vars,'vars');
[L.dom,R.dom] = split_LR(dom,'dom');
[L.depmat,R.depmat] = split_LR(depmat,'depmat');
L.vars = reshape(L.vars,1,[]);      R.vars = reshape(R.vars,1,[]);
% The function spaces must be specified, so that the primary variables can
% be combined by name below
side_names = {'left','right'};
for k=1:2
    if k==1,    S = L;
    else,       S = R;
    end
    Nk = numel(S.vars);
    if ~iscellstr(S.vars) || numel(unique(S.vars))~=Nk
        error("Variable names of the %s operator must be specified as cellstr of distinct names.",side_names{k})
    elseif ~isequal(size(S.dom),[Nk,2])
        error("Domains of the %s operator must be specified as N x 2 array, for N variables.",side_names{k})
    elseif ~isequal(size(S.depmat),[numel(S.dims),Nk])
        error("Dependency array of the %s operator must be specified as p x N array, for p states and N variables.",side_names{k})
    end
end
% Primary variables: those of ZopL, then the new ones of ZopR
[vars_out,idcs_first,idcs_all] = unique([L.vars,R.vars],'stable');
dom_all = [L.dom;R.dom];
out = struct('vars',{reshape(vars_out,1,[])},'dom',dom_all(idcs_first,:));
if any(any(dom_all~=out.dom(idcs_all,:)))
    error("Variables appearing in both the left and right operator must be defined on the same domain.")
end
N = numel(out.vars);
% Map all states to the same function space for the inner product, so that
% Pmat(s) acts on the output of both ZopL and ZopR
is_used = false(1,N);       % variables on which some state depends
is_used(idcs_all([any(L.depmat,1),any(R.depmat,1)])) = true;
L.depmat_out = repmat(is_used,numel(L.dims),1);
R.depmat_out = repmat(is_used,numel(R.dims),1);

% Set default monomial degrees of 0 in the independent variables
if nargin<=5
    degs = 0;
end
if nargin<=6
    opts = struct();    % defaults of monopmat
end

% Pmat is a general matrix by default, symmetric if opts.sym, and positive
% semidefinite (so symmetric) if opts.pos
if ~isstruct(opts)
    error("Options must be specified as a 'struct'.")
end
check_opts(opts,'piesos_quadopvar');
is_pos = isfield(opts,'pos') && opts.pos;
is_sym = is_pos || (isfield(opts,'sym') && opts.sym);


% Check if a psatz multiplier is desired
% --> only supported for positive operator
psatz_is_log = isfield(opts,'psatz') && is_pos && islogical(opts.psatz);
if ~isfield(opts,'psatz') || ~is_pos
    % A single term with multiplier g=1 (find(0) would give no term at all)
    psatz = true;
    npsatz = ones(1,N);
elseif isnumeric(opts.psatz)
    if any(abs(round(opts.psatz))~=opts.psatz)
        error("Maximal order of psatz terms must be specified as nonnegative integers.")
    end
    if isscalar(opts.psatz)
        % Assume same order of psatz along each variable on which a state depends
        npsatz = ones(1,N);
        npsatz(is_used) = opts.psatz+1;
        psatz = true([npsatz,1,1]);
    elseif numel(opts.psatz)==N
        % A separate order of psatz for each variable
        npsatz = reshape(opts.psatz,1,[])+1;
        psatz = true([npsatz,1,1]);
    else
        error("'opts.psatz' must be specified as 1 x N array of integers or ND logical array.")
    end
    opts = rmfield(opts,'psatz');
% Logical array: the psatz terms used, one dimension per variable
elseif islogical(opts.psatz)
    if ~any(opts.psatz,'all')
        error("At least one of the elements of 'opts.psatz' must be true.")
    end
    if N==1 && isvector(opts.psatz)
        opts.psatz = opts.psatz(:);
    end
    npsatz = size(opts.psatz,1:N);
    if numel(opts.psatz)~=prod(npsatz)
        error("'opts.psatz' must be specified as 1 x N array of integers or ND logical array.")
    end
    psatz = opts.psatz;
    opts = rmfield(opts,'psatz');
else
    error("'opts.psatz' must be specified as 1 x N array of integers or ND logical array.")
end
% A multiplier in a variable on which no state depends has no meaning
if any(npsatz(~is_used)>1)
    error("Psatz multipliers cannot be used in variables on which no state depends.")
end
% Options of ZopL, ZopR: those in opts.L, opts.R, the same for both if
% only one is specified, and those in opts itself (exclude, sep, otype)
% for each operator for which opts.L or opts.R does not specify them
if isfield(opts,'L') && ~isfield(opts,'R')
    opts.R = opts.L;
elseif isfield(opts,'R') && ~isfield(opts,'L')
    opts.L = opts.R;
elseif ~isfield(opts,'L')
    opts.L = struct();      opts.R = struct();
end
if ~isstruct(opts.L) || ~isstruct(opts.R)
    error("Options 'L' and 'R' for the operators ZopL and ZopR must be specified as 'struct'.")
end
L.opts = opts.L;    R.opts = opts.R;
for fname = {'exclude','sep','otype'}
    if isfield(opts,fname{1}) && ~isfield(L.opts,fname{1})
        L.opts.(fname{1}) = opts.(fname{1});
    end
    if isfield(opts,fname{1}) && ~isfield(R.opts,fname{1})
        R.opts.(fname{1}) = opts.(fname{1});
    end
end

% Split the degrees of the primary (output) and dummy (input) variables, so
% that Zop has kernels in the dummy variables only, and the degrees of the
% primary variables are used to parameterize Pmat as an SOS matrix below.
% Use the same degrees for both, unless L or R is a field of degs
L.degs = degs;      R.degs = degs;
if isstruct(degs) && (isfield(degs,'L') || isfield(degs,'R'))
    if ~isfield(degs,'L')
        degs.L = degs.R;
    elseif ~isfield(degs,'R')
        degs.R = degs.L;
    end
    L.degs = degs.L;    R.degs = degs.R;
end
% Degrees of the operators ZopL, ZopR: those of the input, with degs.out=0.
for k=1:2
    if k==1,    S = L;
    else,       S = R;
    end
    if isnumeric(S.degs) || iscell(S.degs)
        S.degs_Zop = struct('in',{S.degs},'out',0);
    else
        S.degs_Zop = S.degs;
        S.degs_Zop.out = 0;
    end
    if k==1,    L = S;
    else,       R = S;
    end
end

% Declare a basis of operators on this domain
% Both map to the space of the primary variables. If the inputs of the two
% coincide, ZopR=ZopL
Zop_of = @(S) monopmat(S.dims,struct('in',{S.vars},'out',{out.vars}),...
                    struct('in',S.dom,'out',out.dom),...
                    struct('in',S.depmat,'out',S.depmat_out),S.degs_Zop,S.opts);
ZopL = Zop_of(L);
if isequal({L.dims,L.vars,L.dom,L.depmat,L.degs_Zop,L.opts},...
            {R.dims,R.vars,R.dom,R.depmat,R.degs_Zop,R.opts})
    ZopR = ZopL;
else
    ZopR = Zop_of(R);
end
% For each block row of ZopL and ZopR, the block (state) on which
% it acts, and its number of copies of the primary variables (1 for a
% pointwise output). With otype=3, two block rows act on the same block
[ridcs,cidcs] = find(~cellfun(@isempty,ZopL.ops));
L.blk_idcs(ridcs,1) = cidcs;
L.nfac = max([full(ZopL.depmat.out),ones(numel(ridcs),1)],[],2);
[ridcs,cidcs] = find(~cellfun(@isempty,ZopR.ops));
R.blk_idcs(ridcs,1) = cidcs;
R.nfac = max([full(ZopR.depmat.out),ones(numel(ridcs),1)],[],2);
% A tensor-product block row (otype(1)=2) maps to several copies of the
% primary variables, which are then the variables of Pmat, Zs and psatz.
% Taken from the block rows, not opts.otype: for otype 3, a block of one
% factor is a single row of otype 1, whose Pmat is in si, not si_1
is_tens = any(cellfun(@(K) ~isempty(K) && K.otype(1)==2,[ZopL.ops(:);ZopR.ops(:)]));
if is_tens && psatz_is_log
    error("This case is currently not supported: psatz must be specified as numeric (1 x N array of orders) for otype(1)=2.")
end

% Regions of integration (opts.intmat): rows for the ncop(i) copies of each
% output variable, the most of any block row; a Pmat for each combination
ncop = full(max([ZopL.depmat.out;ZopR.depmat.out],[],1));
intmat = 'full';
if isfield(opts,'intmat')
    intmat = opts.intmat;
end
if ~iscell(intmat)
    intmat = repmat({intmat},1,numel(ncop));
elseif numel(intmat)~=N
    error("Option 'intmat' must be a 1 x N cell, for the N primary variables.")
elseif ~all(cellfun(@(A) (ischar(A) || isstring(A)) && strcmp(A,'full'),intmat(~is_used)))
    error("Regions of integration cannot be specified for a variable on which no state depends.")
else
    % Given per primary variable, as psatz and degs.out: reorder to the
    % output variables of Zop, sorted by name and only those used
    [~,loc_Z] = ismember(ZopL.vars.out,out.vars);
    intmat = intmat(loc_Z);
end
for i=1:numel(ncop)
    intmat{i} = get_intmat(intmat{i},ncop(i));
end
nreg = [cellfun(@(A) size(A,1),intmat),1];
is_full = all(cellfun(@(A,d) isequal(A,repmat([0,d+1],1,d)),intmat,num2cell(ncop)));

% The primary variables of Pmat: the output variables, or for otype(1)=2
% their copies, with the degrees of the monomials of Zs in each of them.
% Declared after ZopL and ZopR, which define the block rows
if is_tens
    vars_P = getvars_tensor(out,L,R,npsatz);
else
    vars_P = getvars_pointwise(out,L,R,npsatz,psatz);
end

% Declare psatz multipliers g_j(s_i) = (s_i-a_i)^j * (b_i-s_i)^(q_i-j),
% nonnegative on dom, in each primary variable s_i
Np = numel(vars_P.names);
pvars = polynomial(reshape(vars_P.names,[],1));
doms = vars_P.dom;
g_psatz = cell(1,Np);
for i=1:Np
    for j=0:vars_P.npsatz(i)-1
        % Store at j+1 (j=0 is not a valid index); the exponents of (s-a)
        % and (b-s) must add up to the order qi = npsatz(i)-1
        g_psatz{i}{j+1} = (pvars(i)-doms(i,1))^(j)*...
                            (doms(i,2)-pvars(i))^(vars_P.npsatz(i)-1-j);
    end
end

% Declare for each block row i of ZopL, ZopR a basis ZsL{i}, ZsR{i} of
% monomials in the primary variables, with exponents degs_out{i} of
% vars_P.L, vars_P.R (variable 1 slowest). In si,
% exponents above max(degs_out)-ceil(qi/2) are cut, for psatz order qi, so
% that each term g_k*Zs'*Qk*Zs has degree <= 2*max(degs_out).
red_deg = ceil((vars_P.npsatz-1)/2);
ZsL = declare_Zs(vars_P.L.degs_out,vars_P.L.dep,vars_P.names,red_deg);
ZsR = declare_Zs(vars_P.R.degs_out,vars_P.R.dep,vars_P.names,red_deg);
PdimL = ZopL.dims.out(:);     % number of rows in each block row of ZopL
PdimR = ZopR.dims.out(:);     % number of rows in each block row of ZopR

% Check symmetric option.
% Pmat(s)=ZsL'*Q*ZsR is symmetric, for symmetric Q, only if the blocks
% match: same number of block rows and columns of the same size, and the
% same monomials in the primary variables. This does not require ZopL=ZopR
if is_sym
    if numel(PdimL)~=numel(PdimR) || any(PdimL~=PdimR)
        error("Output dimensions of ZopL and ZopR must match to declare Pmat as symmetric or positive semidefinite.")
    elseif ~isequal(vars_P.L.degs_out,vars_P.R.degs_out)
        error("Degrees of the primary variables (degs.L.out and degs.R.out) must match to declare Pmat as symmetric or positive semidefinite.")
    end
end
% Option of sosquadvar: Q general, symmetric or positive semidefinite
sos_opt = {};
if is_pos
    sos_opt = {'pos'};
elseif is_sym
    sos_opt = {'sym'};
end

% The psatz multipliers of each region, g_reg{r}: for the whole domain,
% the products of the multipliers g_psatz of the primary variables; for
% other regions, the products of the faces of the region (region_psatz)
nP = prod(nreg);
g_reg = cell([nreg,1]);
if is_full
    g_reg{1} = {1};         % no psatz multiplier: g=1 (also covers N=0)
    lidcs = find(vars_P.psatz);         % linear indices of the terms
    if any(vars_P.npsatz>1)
        g_reg{1} = cell(1,numel(lidcs));
        for t=1:numel(lidcs)
            % Determine order of psatz associated with index lidx
            idcs = num2cell(ones(1,Np+1));
            [idcs{:}] = ind2sub([vars_P.npsatz,1],lidcs(t));
            idcs = cell2mat(idcs(1:Np));
            % Form the psatz multiplier in all variables, once for all blocks
            g_reg{1}{t} = 1;
            for i=1:Np
                g_reg{1}{t} = g_reg{1}{t}*g_psatz{i}{idcs(i)};
            end
        end
    end
else
    [~,loc] = ismember(ZopL.vars.out,out.vars);     % npsatz follows out.vars
    rows = cell(1,numel(nreg));
    for r=1:nP
        [rows{:}] = ind2sub(nreg,r);
        g_reg{r} = {1};
        for i=find(ncop>0)
            cps = get_tensorvars(ZopL.vars.out{i},ncop(i),{'full','pvar'});
            g_i = region_psatz(intmat{i}(rows{i},:),cps,ZopL.dom.out(i,:),npsatz(loc(i))-1);
            [idcs1,idcs2] = ndgrid(1:numel(g_reg{r}),1:numel(g_i));
            g_reg{r} = arrayfun(@(k1,k2) g_reg{r}{k1}*g_i{k2},idcs1(:).',idcs2(:).','UniformOutput',false);
        end
    end
end

% For each psatz term, declare an SOS matrix acting on the operator basis,
% so that Pop = ZopL'*(sum_k g_k*Pmat_k)*ZopR, for Pmat_k = ZsL'*Qk*ZsR, Qk>=0.
% The matrix Pmat_k is a p x q cell, a block for each pair of rows in ZopL.ops
% and ZopR.ops. One such Pmat for each region, Pmat{r}
Pmat = cell([nreg,1]);
for r=1:nP
    g_terms = g_reg{r};
    if isequal(g_terms,{1})
        % No psatz multiplier: Q general, symmetric or positive semidefinite
        [prog,Pmat{r}] = sosquadvar(prog,ZsL,ZsR,PdimL,PdimR,sos_opt{:});
        continue
    end
    Pmat_terms = cell(1,numel(g_terms));
    for t=1:numel(g_terms)
        % Declare an SOS matrix
        [prog,Pmat_terms{t}] = sosquadvar(prog,ZsL,ZsR,PdimL,PdimR,'pos');
    end
    % Add the terms to the full nonnegative matrix variable, for each block
    % in a single product [g_1*I,...,g_T*I]*[Pmat_1;...;Pmat_T].
    q = numel(PdimL);       % psatz requires opts.pos, so p=q
    Pmat{r} = cell(q,q);
    for l=1:q^2
        row_l = mod(l-1,q)+1;       % block row of Pmat{r}{l}
        Pstack = cellfun(@(P) P{l},Pmat_terms,'UniformOutput',false);
        Gstack = cellfun(@(g) g*speye(PdimL(row_l)),g_terms,'UniformOutput',false);
        Pmat{r}{l} = horzcat(Gstack{:})*vertcat(Pstack{:});
    end
end

% A single Pmat for the whole domain
if is_full
    Pmat = Pmat{1};
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [XL,XR] = split_LR(X,name)
% [XL,XR] = SPLIT_LR(X,NAME) extracts the variables, or dimensions,
% or domains, or depmats for the left and right operator from the array or
% struct X, with fields L (or out) and/or R (or in). If X is an array, or
% if only one of the fields is specified, the same value is used for both.
% NAME is used in error messages.
%
% INPUTS
% - X:      array, or 'struct' with fields L (or out) and/or R (or in);
% - name:   char, the name of the input, for error messages;
%
% OUTPUTS
% - XL, XR: the values for the left and right operator.

% An array is used for both
if ~isstruct(X)
    XL = X;     XR = X;
    return
end
% Each side given at most once, by either name
isL = isfield(X,{'L','out'});       isR = isfield(X,{'R','in'});
if ~any(isL) && ~any(isR)
    error("Input '%s' must be specified as an array, or as a struct with fields L/out and/or R/in.",name)
elseif all(isL) || all(isR)
    error("Input '%s' must specify each of the left (L/out) and right (R/in) values only once.",name)
end
fnamesL = {'L','out'};              fnamesR = {'R','in'};
if any(isL)
    XL = X.(fnamesL{isL});
end
if any(isR)
    XR = X.(fnamesR{isR});
end
% A side not given takes the value of the other
if ~any(isL)
    XL = XR;
elseif ~any(isR)
    XR = XL;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Zs = declare_Zs(degs_out,is_dep_out,varnames,red_deg)
% ZS = DECLARE_ZS(DEGS_OUT,IS_DEP_OUT,VARNAMES,RED_DEG) declares for each
% block row i of an operator basis a basis Zs{i} of monomials in the
% primary variables, with exponents degs_out{i} as returned by split_degs
% or copy_degs. The degree in each variable is cut by the entry of RED_DEG
% for that variable, but to at least the minimal exponent.
%
% INPUTS
% - degs_out:   1 x q cell; degs_out{i} is a cell of column vectors of
%               exponents, one per variable on which block row i depends;
% - is_dep_out: q x P logical array, the primary variables of each block
%               row;
% - varnames:   1 x P cellstr, the names of the primary variables;
% - red_deg:    1 x P array, the degree reduction in each variable;
%
% OUTPUTS
% - Zs:         q x 1 cell; Zs{i} is the nZ_i x 1 'polynomial' basis of
%               block row i, or 1 if it depends on no variable.

q = numel(degs_out);
Zs = cell(q,1);
for i=1:q
    idcs_i = find(is_dep_out(i,:));
    degmat_i = zeros(1,0);
    nZ_i = 1;
    for k=1:numel(idcs_i)
        % Cut the degree in the kth variable, but keep at least min(e_k)
        e_k = degs_out{i}{k};
        e_k = e_k(e_k<=max(max(e_k)-red_deg(idcs_i(k)),min(e_k)));
        degmat_i = [repelem(degmat_i,numel(e_k),1),repmat(e_k,nZ_i,1)];
        nZ_i = nZ_i*numel(e_k);
    end
    if isempty(idcs_i)
        Zs{i} = 1;      % constant basis, as accepted by sosquadvar
    else
        Zs{i} = polynomial(speye(nZ_i),degmat_i,varnames(idcs_i),[nZ_i,1]);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function vars_P = getvars_pointwise(out,L,R,npsatz,psatz)
% VARS_P = GETVARS_POINTWISE(OUT,L,R,NPSATZ,PSATZ) declares the primary
% variables of Pmat for a pointwise output of Zop (otype(1)=1): the output
% variables themselves, with the degrees of the monomials of Zs in each.
%
% INPUTS
% - out:    struct with fields 'vars' (1 x N cellstr, the names of the
%           output variables) and 'dom' (N x 2 array, their domains);
% - L, R:   structs for ZopL and ZopR, with fields
%     depmat_out: p x N logical array, on which output variables each
%                 block depends;
%     degs:       monomial degrees, as specified in piesos_quadopvar;
% - npsatz: 1 x N array, the order+1 of the psatz multipliers in each
%           output variable;
% - psatz:  N-dimensional logical array, true for the multipliers to use;
%
% OUTPUTS
% - vars_P:     struct as returned by getvars_tensor, with a primary
%               variable for each output variable (no copies);

degs_outL = split_degs(L.degs,L.depmat_out);
degs_outR = split_degs(R.degs,R.depmat_out);
vars_P = struct('names',{out.vars},'dom',out.dom,'npsatz',npsatz, ...
    'psatz',psatz, ...
    'L',struct('dep',L.depmat_out,'degs_out',{degs_outL}), ...
    'R',struct('dep',R.depmat_out,'degs_out',{degs_outR}));

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function vars_P = getvars_tensor(out,L,R,npsatz)
% VARS_P = GETVARS_TENSOR(OUT,L,R,NPSATZ) declares the primary variables
% of Pmat for a tensor-product output of Zop (otype(1)=2): the copies si_j,
% j=1,...,d, of each output variable si, with d the largest number of
% copies of a block row. A block row with di copies depends on the copies
% 1,...,di of the output variables on which it depends.
%
% INPUTS
% - out:    struct with fields 'vars' (1 x N cellstr, the names of the
%           output variables) and 'dom' (N x 2 array, their domains);
% - L, R:   structs for ZopL and ZopR, with fields
%     blk_idcs:   r x 1 array, the block (state) on which each block row
%                 of Zop acts;
%     nfac:       r x 1 array, the number of copies of each block row: the
%                 number of factors for a tensor-product output, and 1 for
%                 a pointwise output, or a block without spatial variables;
%     depmat_out: p x N logical array, on which output variables each
%                 block depends;
%     degs:       monomial degrees, as specified in piesos_quadopvar, with
%                 the degrees degs.out of Zs in each copy, see copy_degs;
% - npsatz: 1 x N array, the order+1 of the psatz multipliers in each
%           output variable, used for each of its copies;
%
% OUTPUTS
% - vars_P:     struct with fields
%     names:    1 x P cellstr, the names si_j of the P=N*d copies, ordered
%               by variable and then by copy;
%     dom:      P x 2 array, the domains of the copies;
%     npsatz:   1 x P array, the order+1 of the multiplier in each copy;
%     psatz:    P-dimensional logical array, true for all multipliers;
%     L, R:     structs for ZopL and ZopR, with fields
%         dep:      r x P logical array, on which copies each block row
%                   depends;
%         degs_out: 1 x r cell, with degs_out{i} a 1 x nnz(dep(i,:)) cell
%                   of column vectors with the exponents of Zs{i} in each
%                   copy on which block row i depends, in the order of
%                   names, as for split_degs;

N = numel(out.vars);
sides = {L,R};
% The number of copies of each block row, from ZopL and ZopR
nfac = {L.nfac,R.nfac};
d = max([nfac{1};nfac{2}]);
% The copies, ordered by variable and then by copy
P = N*d;
names = get_tensorvars(out.vars,d,'full');      dom = zeros(P,2);
copies = struct('var',zeros(1,P),'copy',zeros(1,P));
for i=1:N
    for j=1:d
        p = (i-1)*d+j;
        dom(p,:) = out.dom(i,:);
        copies.var(p) = i;      copies.copy(p) = j;
    end
end
npsatz_copies = npsatz(copies.var);
vars_P = struct('names',{names},'dom',dom,'npsatz',npsatz_copies, ...
    'psatz',true([npsatz_copies,1,1]));
for k=1:2
    S = sides{k};
    % A row per block row of Zop, with the dependence of its block
    dep = S.depmat_out(S.blk_idcs,copies.var) & (copies.copy<=nfac{k});
    degs_out = copy_degs(S.degs,S.depmat_out,nfac{k},dep,copies,S.blk_idcs);
    side = struct('dep',dep,'degs_out',{degs_out});
    if k==1
        vars_P.L = side;
    else
        vars_P.R = side;
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function degs_out = copy_degs(degs,depmat_out,nfac,dep,copies,blk_idcs)
% DEGS_OUT = COPY_DEGS(DEGS,DEPMAT_OUT,NFAC,DEP,COPIES,BLK_IDCS) splits
% the monomial degrees DEGS into the degrees of the copies of the output
% variables for each block row of Zop, for the basis Zs of Pmat
% (otype(1)=2, see getvars_tensor).
%
% INPUTS
% - degs:       monomial degrees: numeric, cell, or struct with fields 'in'
%               and 'out'. The degrees degs.out of block i (the same for
%               all blocks, unless a cell with an element for each) can be
%               - an array: a scalar for all copies of all variables, or an
%                 array with a column for each of the Mi output variables
%                 of the block (or for all M), and one row, for all copies
%                 of a variable, or a row for each copy (the first nfac(i)
%                 rows are used if there are more);
%               - a cell of the same shape, of column vectors of exponents;
% - depmat_out: q x M logical array, on which output variables each block
%               depends;
% - nfac:       r x 1 array, the number of copies of each block row;
% - dep:        r x P logical array, on which copies each block row
%               depends;
% - copies:     struct with fields 'var' and 'copy', 1 x P arrays with the
%               output variable and the copy of each of the P copies;
% - blk_idcs:   r x 1 array, the block (1,...,q) of each block row;
%
% OUTPUTS
% - degs_out:   1 x r cell, with degs_out{i} a 1 x nnz(dep(i,:)) cell of
%               column vectors with the exponents for each copy of block
%               row i, in the order of the copies;

[q,M] = size(depmat_out);
if isnumeric(degs) || iscell(degs)
    degs = struct('in',{degs},'out',{degs});
elseif ~isstruct(degs) || ~isfield(degs,'in') || ~isfield(degs,'out')
    error("Options 'degs' must be specified as a struct with fields 'in' and 'out'.")
end
degs_raw = degs.out;
if isnumeric(degs_raw)
    degs_raw = repmat({degs_raw},1,q);
elseif isscalar(degs_raw)
    degs_raw = repmat(degs_raw,1,q);
elseif numel(degs_raw)~=q
    error("Degrees must be specified as integer array or 1 x q cell for q inputs.")
end
% Use the degrees and dependence of the block of each block row
degs_raw = degs_raw(blk_idcs);      depmat_out = depmat_out(blk_idcs,:);
r = numel(blk_idcs);
degs_out = cell(1,r);
for i=1:r
    % Degrees in the Mi output variables of block row i (none if Mi = 0)
    Mi = nnz(depmat_out(i,:));
    if Mi==0
        degs_out{i} = cell(1,0);
        continue
    end
    % A maximal degree a gives the exponents 0..a; a cell lists them
    degval = degs_raw{i};
    if isnumeric(degval)
        if isscalar(degval)
            degval = degval*ones(1,Mi);
        elseif size(degval,2)==M
            degval = degval(:,depmat_out(i,:));
        end
        if any(round(abs(degval(:)))~=degval(:))
            error("Monomial degrees must be specified as nonnegative integers.")
        end
        degval = cellfun(@(a)(0:a)',num2cell(degval),'UniformOutput',false);
    elseif size(degval,2)==M
        degval = degval(:,depmat_out(i,:));
    end
    % A column for each variable, a row for each copy
    if size(degval,2)~=Mi
        error("Degrees of the output variables must be specified as a 1 x M array or cell, for M variables on which the state depends.")
    elseif size(degval,1)==1
        degval = repmat(degval,nfac(i),1);
    elseif size(degval,1)>=nfac(i)
        degval = degval(1:nfac(i),:);
    else
        error("Row dimension of the degrees of the output variables must be 1, or at least the number of factors of the state.")
    end
    % The exponents of each copy on which block row i depends, in the order
    % of the copies. The column of a variable among those of the block is
    % col(variable)
    col = cumsum(depmat_out(i,:));
    idcs = find(dep(i,:));
    degs_out{i} = cell(1,numel(idcs));
    for t=1:numel(idcs)
        degs_out{i}{t} = reshape(degval{copies.copy(idcs(t)),col(copies.var(idcs(t)))},[],1);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function degs_out = split_degs(degs,is_dep_out)
% DEGS_OUT = SPLIT_DEGS(DEGS,IS_DEP_OUT) extracts from the monomial
% degrees DEGS the degrees of the output (primary) variables.
%
% INPUTS
% - degs:       monomial degrees: numeric, cell,
%               or struct with fields 'in' and 'out'. If numeric or cell,
%               the same degrees are used for the input and output;
% - is_dep_out: q x M logical array, specifying for each of the q blocks
%               on which of the M output variables it depends;
%
% OUTPUTS
% - degs_out:   1 x q cell, with degs_out{i} a 1 x Mi cell of column
%               vectors specifying the exponents of the monomials in each
%               of the Mi output variables on which block i depends
%               (maximal degrees d are expanded to (0:d)');
%
% NOTES
% - degs.out is processed for a single factor, as the primary monomial
%   appears only in the first factor of the tensor-PI operator
%   (opts.otype(1)==1);

[q,M] = size(is_dep_out);
if isnumeric(degs) || iscell(degs)
    % Same degrees for input and output
    degs = struct('in',{degs},'out',{degs});
elseif ~isstruct(degs) || ~isfield(degs,'in') || ~isfield(degs,'out')
    error("Options 'degs' must be specified as a struct with fields 'in' and 'out'.")
end

% Extract the degrees of the output variables for each block
degs_out = degs.out;
if isnumeric(degs_out)
    degs_out = repmat({degs_out},1,q);
elseif isscalar(degs_out)
    degs_out = repmat(degs_out,1,q);
elseif numel(degs_out)~=q
    error("Degrees must be specified as integer array or 1 x q cell for q inputs.")
end
for i=1:q
    Mi = nnz(is_dep_out(i,:));
    degval = degs_out{i};
    if Mi==0
        % The block does not depend on any variable
        degs_out{i} = cell(1,0);
        continue
    end
    if isnumeric(degval)
        % Maximal degrees are specified --> include all lower degrees
        if isscalar(degval)
            degval = degval*ones(1,Mi);
        elseif size(degval,2)==M
            degval = degval(:,is_dep_out(i,:));
        end
        if any(round(abs(degval(:)))~=degval(:))
            error("Monomial degrees must be specified as nonnegative integers.")
        end
        degval = cellfun(@(a)(0:a)',num2cell(degval),'UniformOutput',false);
    elseif size(degval,2)==M
        % Retain only monomials in relevant variables
        degval = degval(:,is_dep_out(i,:));
    end
    if ~isequal(size(degval),[1,Mi])
        error("Degrees of the output variables must be specified as a 1 x M array or cell, for M variables on which the state depends.")
    end
    degs_out{i} = cellfun(@(a)a(:),degval,'UniformOutput',false);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function intmat = get_intmat(spec,d)
% INTMAT = GET_INTMAT(SPEC,D) returns the regions SPEC of the copies
% s_1,...,s_D of a variable as the intmat of an 'intop', a row per region,
% checking a given array, or building a preset (see opts.intmat).
%
% INPUTS
% - spec:   q x 2*D array of integers, with the limits of s_j among 0 (a),
%           D+1 (b) and 1,...,j-1 (s_1,...,s_(j-1)), or a preset 'full',
%           'steps' or 'all';
% - d:      nonnegative integer, the number of copies D;
%
% OUTPUTS
% - intmat: q x 2*D array, a row per region;

if isnumeric(spec)
    lim_max = repelem(0:d-1,2);     % the last copy in the limits of each copy
    if isempty(spec) || size(spec,2)~=2*d || any(spec~=round(spec) | spec<0 | spec>d+1 | (spec>lim_max & spec<=d),'all')
        error("Regions of %d copies must be specified as q x %d array, with the limits of copy j among 0 (a), %d (b) and the copies 1,...,j-1.",d,2*d,d+1)
    end
    intmat = spec;
    return
elseif ~(ischar(spec) || isstring(spec)) || ~ismember(spec,{'full','steps','all'})
    error("Presets of option 'intmat' are 'full', 'steps' and 'all'.")
end
intmat = repmat([0,d+1],1,d);
if strcmp(spec,'full') || d<=1
    return
end
% Insert the copies one at a time. An ordering of s_1,...,s_(j-1) is stored
% as ords{k}, the first copy of each group of tied copies, in increasing
% order, with ntie(k) ties; s_j enters between two groups (or at an end),
% or joins a group (a tie, the Dirac delta [o,o])
ords = {1};     ntie = 0;
for j=2:d
    rows = zeros(0,2*d);    ords_new = cell(1,0);   ntie_new = zeros(1,0);
    for k=1:numel(ords)
        o = ords{k};        lims = [0,o,d+1];
        for p=1:numel(o)+1
            rows(end+1,:) = intmat(k,:);                                    %#ok<AGROW>
            rows(end,2*j-1:2*j) = lims([p,p+1]);
            ords_new{end+1} = [o(1:p-1),j,o(p:end)];                        %#ok<AGROW>
            ntie_new(end+1) = ntie(k);                                      %#ok<AGROW>
        end
        for p=1:numel(o)
            rows(end+1,:) = intmat(k,:);                                    %#ok<AGROW>
            rows(end,2*j-1:2*j) = o([p,p]);
            ords_new{end+1} = o;                                            %#ok<AGROW>
            ntie_new(end+1) = ntie(k)+1;                                    %#ok<AGROW>
        end
    end
    intmat = rows;      ords = ords_new;        ntie = ntie_new;
end
% Strict orderings only, or all orderings except all copies tied (d-1 ties)
if strcmp(spec,'steps')
    intmat = intmat(ntie==0,:);
else
    intmat = intmat(ntie<d-1,:);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function g = region_psatz(introw,var_cps,dom,q)
% G = REGION_PSATZ(INTROW,VAR_CPS,DOM,Q) returns the psatz multipliers of a
% region of the copies of a variable: for each group of k copies related by
% the limits of the region, the products of total order k*q of the faces of
% the group, and the products of these over the groups. The faces are the
% differences y-x of values x<y with no value between them, among a, b and
% the copies (a tied copy taking the value of its target). A copy over
% [a,b] is a group by itself, with faces s-a and b-s, as on the whole
% domain; the lower triangle s_2<=s_1 is one group, with faces s_2-a,
% s_1-s_2 and b-s_1.
%
% INPUTS
% - introw:  1 x 2*d array, a region of the d copies, a row of opts.intmat;
% - var_cps: d x 1 'polynomial', the copies s_1,...,s_d;
% - dom:     1 x 2 array [a,b], the domain of the variable;
% - q:       nonnegative integer, the order of psatz (in each copy);
%
% OUTPUTS
% - g:       1 x T cell of the multipliers, nonnegative on the region ({1}
%            for q=0);

g = {1};
if q==0
    return
end
% Values (nodes): a (1), the copies s_1,...,s_d (2,...,d+1), b (d+2). A
% tied copy takes the node of its target, an earlier copy
d = numel(var_cps);     n = d+2;
L = introw(1:2:end)+1;     U = introw(2:2:end)+1;
node = 1:n;
for j=find(L==U)
    node(j+1) = node(L(j));
end
% lt(x,y): x<y on the region, from a<s_j<b and the limits of the copies
lt = false(n);      lt(1,2:n) = true;       lt(1:n-1,n) = true;
for j=find(L~=U)
    lt(node(L(j)),j+1) = true;      lt(j+1,node(U(j))) = true;
end
for k=1:n
    lt = lt | (double(lt)*double(lt))>0;
end
is_used = false(1,n);   is_used([1,node(2:n-1),n]) = true;
lt(~is_used,:) = false;     lt(:,~is_used) = false;     lt(logical(eye(n))) = false;
% Faces: x<y with no value between them
[xs,ys] = find(lt & ~((double(lt)*double(lt))>0));
vals = [{dom(1)},arrayfun(@(j) var_cps(j),1:d,'UniformOutput',false),{dom(2)}];
% Groups of copies, connected by faces between two copies
grp = 1:n;
is_cc = xs>1 & xs<n & ys>1 & ys<n;
for k=1:n
    for f=find(is_cc).'
        grp([xs(f),ys(f)]) = min(grp([xs(f),ys(f)]));
    end
end
free = find(is_used(2:n-1))+1;
for c=unique(grp(free))
    C = free(grp(free)==c);
    faces = find(ismember(xs,C) | ismember(ys,C));
    % Exponents of the faces with total order k*q, for k copies in the group
    m = numel(C)*q;     E = zeros(1,0);
    for f=1:numel(faces)
        E = [repelem(E,m+1,1),repmat((0:m)',size(E,1),1)];
        E = E(sum(E,2)<=m,:);
    end
    E = E(sum(E,2)==m,:);
    g_c = cell(1,size(E,1));
    for e=1:size(E,1)
        g_c{e} = 1;
        for f=find(E(e,:)>0)
            g_c{e} = g_c{e}*(vals{ys(faces(f))}-vals{xs(faces(f))})^E(e,f);
        end
    end
    [idcs1,idcs2] = ndgrid(1:numel(g),1:numel(g_c));
    g = arrayfun(@(k1,k2) g{k1}*g_c{k2},idcs1(:).',idcs2(:).','UniformOutput',false);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
