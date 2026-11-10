function [prog,Pmat,ZopL,ZopR] = piesos_quadopvar(prog,dims,vars,dom,depmat,degs,opts)
% [PROG,PMAT,ZOPL,ZOPR] = PIESOS_QUADOPVAR(PROG,DIMS,VARS,DOM,DEPMAT,
% DEGS,OPTS) declares an operator variable of the form
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
%           operator variable (see `piesos_program`);
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
%           variables in each operator block. Each element degs.in{i} can
%           be a di x (M or N) cell, where di depends on depmat.L/R(i,:),
%           specifying the degrees of the monomial basis in each of the
%           factors of each of the input variables in block i. Element
%           in{i} can also be a di x (M or N) array, specifying the
%           maximal monomial degree in each of the factors of each
%           variable. If in{i} is scalar, the same degree is used for all
%           variables;
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
%           ZopL and ZopR. Defaults to 0;
% - opts:   (optional) struct with fields
%     pos:  set to true to declare Pmat as a positive semidefinite (SOS)
%           variable. Requires output dimensions of ZopL and ZopR to
%           match (note that output dimensions depend on monomial
%           degrees), and the degrees degs.L.out and degs.R.out of the
%           primary variables to match. Defaults to false;
%     sym:  set to true to declare Pmat as a symmetric variable
%           (Pmat(s)'=Pmat(s)). Requires output dimensions of ZopL and
%           ZopR to match (note that output dimensions depend on
%           monomial degrees), and degs.L.out and degs.R.out to match.
%           Defaults to false;
%   psatz:  (only if opts.pos==true) 1 x P array of integers, specifying
%           the maximal order of the psatz multiplier in each of the P
%           primary variables. Here, psatz(i)=qi implies inclusion of
%           all multipliers
%               gj = (bi-si)^{qi-j}*(si-ai)^{j},    j in {0,...,qi}
%           This means a total number of prod(psatz+1) terms. Can also
%           be set to a q1 x ... x qP logical array, where now the
%           dimensions of psatz determine the maximal order, and only the
%           multipliers associated with true elements of psatz are
%           included. Element k along dimension i of psatz corresponds to
%           gj with j=k-1, so psatz(i)=qi is equivalent to a logical
%           array of size qi+1 along dimension i. Terms are ordered by
%           linear index in psatz. A scalar psatz is used for all
%           variables on which some state depends. For gj of order qi,
%           the degree in si of the basis monomials of Pmat is reduced by
%           ceil(qi/2), so that the degree of Pmat in si does not exceed
%           twice its degree in degs.out;
%    L/R:   (optional) fields specifying separate options for the
%           operators ZopL and ZopR. Each is a struct with fields
%           'exclude', 'sep', and 'otype'. See `monopmat` function for
%           more details on how these options can be specified (the case
%           otype(1)==2 is not supported). If only one of L, R is
%           specified, it is used for both. If neither is, the fields
%           'exclude', 'sep', 'otype' of opts are used for both;
%
% OUTPUTS
% - prog:   'struct' specifying the PIESOS program with the matrix
%           variable Pmat added;
% - Pmat:   p x q cell, with element Pmat{i,j} an mi x nj 'dpvar' object
%           in the primary variables, where mi and nj are the numbers of
%           rows of the ith block row of ZopL.ops and the jth block row
%           of ZopR.ops (p=q if opts.sym or opts.pos). The cell
%           represents the matrix Pmat(s) parameterizing
%           Pop = ZopL'*Pmat*ZopR, with Pmat{i,j} of the form ZsL'*Q*ZsR
%           for bases ZsL, ZsR of monomials in the primary variables, and
%           Q general, symmetric (opts.sym) or positive semidefinite
%           (opts.pos, an SOS matrix). With psatz multipliers g_k,
%           Pmat = sum_k g_k*Pmat_k, where Pmat_k is of the form
%           Zs'*Qk*Zs, with Qk >= 0 and Zs a basis of monomials in the
%           primary variables;
% - ZopL/ZopR:  rtot x (m or n) 'tensopmat' object representing a basis of
%           operators of specified degrees on the declared function space,
%           where m/n = sum(dims.L/R) and where rtot depends on the
%           selected monomial degrees. The kernels of the operators
%           depend only on the dummy variables;
%
% NOTES
% - For dims, vars, dom, depmat and degs, if specified as an array, or if
%   only one of the fields L, R is specified, the same value is used for
%   ZopL and ZopR.
% - The polynomial degree of Pmat depends on the output degree of both
%   ZopL and ZopR, as Pmat is constructed using the `sosquadvar` function.
% - For states depending on different variables, the operators ZopL and
%   ZopR map to a function space in all of the variables appearing in
%   any state of ZopL or ZopR (the primary variables). For example, for
%   state x(s1) in ZopL and y(s1,s2) in ZopR, both map to functions of s1
%   and s2. The primary variables are ordered as vars.L, followed by those
%   of vars.R not in vars.L; degs.L/R.out and psatz refer to this order.

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
% DJ, 10/15/2026: Initial coding;

% Check that the inputs are appropriate
% ZopL acts on the space (dimsL,varsL,domL,depmatL), ZopR on (dimsR,varsR,
% domR,depmatR).
[dimsL,dimsR] = split_LR(dims,'dims');
[varsL,varsR] = split_LR(vars,'vars');
[domL,domR] = split_LR(dom,'dom');
[depmatL,depmatR] = split_LR(depmat,'depmat');
varsL = reshape(varsL,1,[]);        varsR = reshape(varsR,1,[]);
check_space(dimsL,varsL,domL,depmatL,'left');
check_space(dimsR,varsR,domR,depmatR,'right');
% Primary variables: those of ZopL, then the new ones of ZopR
[vars_out,idcs_first,idcs_all] = unique([varsL,varsR],'stable');
vars_out = reshape(vars_out,1,[]);
dom_all = [domL;domR];
dom_out = dom_all(idcs_first,:);
if any(any(dom_all~=dom_out(idcs_all,:)))
    error("Variables appearing in both the left and right operator must be defined on the same domain.")
end
N = numel(vars_out);
% Map all states to the same function space for the inner product, so that
% Pmat(s) acts on the output of both ZopL and ZopR
is_used = false(1,N);       % variables on which some state depends
is_used(idcs_all([any(depmatL,1),any(depmatR,1)])) = true;
depmatL_out = repmat(is_used,numel(dimsL),1);
depmatR_out = repmat(is_used,numel(dimsR),1);

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
is_pos = isfield(opts,'pos') && opts.pos;
is_sym = is_pos || (isfield(opts,'sym') && opts.sym);


% Check if a psatz multiplier is desired
% --> only supported for positive operator
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
% A psatz multiplier only has meaning for an SOS matrix
if any(npsatz>1) && ~is_pos
    error("Psatz multipliers can only be used if Pmat is positive semidefinite (opts.pos=true).")
end
% Options of ZopL, ZopR: those in opts.L, opts.R, the same for both if only
% one is specified, and those in opts itself if neither is
if isfield(opts,'L') || isfield(opts,'R')
    if ~isfield(opts,'L')
        opts.L = opts.R;
    elseif ~isfield(opts,'R')
        opts.R = opts.L;
    end
    optsL = opts.L;     optsR = opts.R;
else
    optsL = opts;       optsR = opts;
end
if ~isstruct(optsL) || ~isstruct(optsR)
    error("Options 'L' and 'R' for the operators ZopL and ZopR must be specified as 'struct'.")
end
if (isfield(optsL,'otype') && optsL.otype(1)==2) || (isfield(optsR,'otype') && optsR.otype(1)==2)
    error("This case is currently not supported.")
end

% Split the degrees of the primary (output) and dummy (input) variables, so
% that Zop has kernels in the dummy variables only, and the degrees of the
% primary variables are used to parameterize Pmat as an SOS matrix below.
% Use the same degrees for both, unless L or R is a field of degs
degsL = degs;       degsR = degs;
if isstruct(degs) && (isfield(degs,'L') || isfield(degs,'R'))
    if ~isfield(degs,'L')
        degs.L = degs.R;
    elseif ~isfield(degs,'R')
        degs.R = degs.L;
    end
    degsL = degs.L;     degsR = degs.R;
end
[degs_outL,degs_newL] = split_degs(degsL,depmatL_out);
[degs_outR,degs_newR] = split_degs(degsR,depmatR_out);

% Declare a basis of operators on this domain
% Both map to the space of the primary variables. If the inputs of the two
% coincide, ZopR=ZopL
ZopL = monopmat(dimsL,struct('in',{varsL},'out',{vars_out}),...
                struct('in',domL,'out',dom_out),...
                struct('in',depmatL,'out',depmatL_out),degs_newL,optsL);
if isequal({dimsL,varsL,domL,depmatL,degs_newL,optsL},{dimsR,varsR,domR,depmatR,degs_newR,optsR})
    ZopR = ZopL;
else
    ZopR = monopmat(dimsR,struct('in',{varsR},'out',{vars_out}),...
                    struct('in',domR,'out',dom_out),...
                    struct('in',depmatR,'out',depmatR_out),degs_newR,optsR);
end

% Declare psatz multipliers g_j(s_i) = (s_i-a_i)^j * (b_i-s_i)^(q_i-j),
% nonnegative on dom
pvars = polynomial(reshape(vars_out,[],1));
doms = dom_out;
g_psatz = cell(1,N);
for i=1:N
    for j=0:npsatz(i)-1
        % Store at j+1 (j=0 is not a valid index); the exponents of (s-a)
        % and (b-s) must add up to the order qi = npsatz(i)-1
        g_psatz{i}{j+1} = (pvars(i)-doms(i,1))^(j)*...
                            (doms(i,2)-pvars(i))^(npsatz(i)-1-j);
    end
end

% Declare for each block row i of ZopL, ZopR a basis ZsL{i}, ZsR{i} of
% monomials in the primary variables, with exponents degs_outL{i},
% degs_outR{i} (variable 1 slowest, as for ZL in sopvar). In si, exponents
% above max(degs_out)-ceil(qi/2) are cut, for psatz order qi, so that each
% term g_k*Zs'*Qk*Zs has degree <= 2*max(degs_out). 
varnames = vars_out;
red_deg = ceil((npsatz-1)/2);
ZsL = declare_Zs(degs_outL,depmatL_out,varnames,red_deg);
ZsR = declare_Zs(degs_outR,depmatR_out,varnames,red_deg);
PdimL = ZopL.dims.out(:);     % number of rows in each block row of ZopL
PdimR = ZopR.dims.out(:);     % number of rows in each block row of ZopR

% Check symmetric option.
% Pmat(s)=ZsL'*Q*ZsR is symmetric, for symmetric Q, only if the blocks
% match: same number of block rows and columns of the same size, and the
% same monomials in the primary variables. This does not require ZopL=ZopR
if is_sym
    if numel(PdimL)~=numel(PdimR) || any(PdimL~=PdimR)
        error("Output dimensions of ZopL and ZopR must match to declare Pmat as symmetric or positive semidefinite.")
    elseif ~isequal(degs_outL,degs_outR)
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

% For each psatz term, declare an SOS matrix acting on the operator basis,
% so that Pop = ZopL'*(sum_k g_k*Pmat_k)*ZopR, for Pmat_k = ZsL'*Qk*ZsR, Qk>=0.
% The matrix Pmat_k is a p x q cell, a block for each pair of rows in ZopL.ops
% and ZopR.ops.
if all(npsatz==1)
    % No psatz multiplier: single term with g=1 (also covers N=0)
    [prog,Pmat] = sosquadvar(prog,ZsL,ZsR,PdimL,PdimR,sos_opt{:});
else
    lidcs = find(psatz);        % linear indices of the terms
    Pmat_terms = cell(1,numel(lidcs));
    g_terms = cell(1,numel(lidcs));
    for t=1:numel(lidcs)
        % Declare an SOS matrix
        [prog,Pmat_terms{t}] = sosquadvar(prog,ZsL,ZsR,PdimL,PdimR,'pos'); 
        % Determine order of psatz associated with index lidx
        idcs = num2cell(ones(1,N+1));
        [idcs{:}] = ind2sub([npsatz,1],lidcs(t));
        idcs = cell2mat(idcs(1:N));
        % Form the psatz multiplier in all variables, once for all blocks
        g_terms{t} = 1;
        for i=1:N
            g_terms{t} = g_terms{t}*g_psatz{i}{idcs(i)};
        end
    end
    % Add the terms to the full nonnegative matrix variable, for each block
    % in a single product [g_1*I,...,g_T*I]*[Pmat_1;...;Pmat_T].
    q = numel(PdimL);       % psatz requires opts.pos, so p=q 
    Pmat = cell(q,q);
    for l=1:numel(Pmat)
        row_l = mod(l-1,q)+1;       % block row of Pmat{l}
        Pstack = cellfun(@(P) P{l},Pmat_terms,'UniformOutput',false);
        Gstack = cellfun(@(g) g*speye(PdimL(row_l)),g_terms,'UniformOutput',false); 
        Pmat{l} = horzcat(Gstack{:})*vertcat(Pstack{:});
    end
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

if ~isstruct(X)
    XL = X;     XR = X;
    return
end
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
if ~any(isL)
    XL = XR;
elseif ~any(isR)
    XR = XL;
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_space(dims,vars,dom,depmat,side)
% CHECK_SPACE(DIMS,VARS,DOM,DEPMAT,SIDE) checks that the function space of
% the left or right operator (SIDE='left' or 'right') is properly specified,
% so that the primary variables can be combined by name below.

p = numel(dims);        N = numel(vars);
if ~iscellstr(vars) || numel(unique(vars))~=N
    error("Variable names of the %s operator must be specified as cellstr of distinct names.",side)
elseif ~isequal(size(dom),[N,2])
    error("Domains of the %s operator must be specified as N x 2 array, for N variables.",side)
elseif ~isequal(size(depmat),[p,N])
    error("Dependency array of the %s operator must be specified as p x N array, for p states and N variables.",side)
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Zs = declare_Zs(degs_out,is_dep_out,varnames,red_deg)
% ZS = DECLARE_ZS(DEGS_OUT,IS_DEP_OUT,VARNAMES,RED_DEG) declares for each
% block row i of an operator basis a basis Zs{i} of monomials in the primary
% variables, with exponents degs_out{i} as returned by split_degs. The
% degree in each variable is cut by the entry of RED_DEG for that variable,
% but to at least the minimal exponent.

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
function [degs_out,degs_new] = split_degs(degs,is_dep_out)
% [DEGS_OUT,DEGS_NEW] = SPLIT_DEGS(DEGS,IS_DEP_OUT) splits the monomial
% degrees DEGS, specified as in 'monopmat', into the degrees of the output
% (primary) variables, and an updated set of degrees where the degrees of
% the output variables are set to 0.
%
% INPUTS
% - degs:       monomial degrees as specified in 'monopmat': numeric, cell,
%               or struct with fields 'in' and 'out'. If numeric or cell,
%               the same degrees are used for the input and output;
% - is_dep_out: q x M logical array, specifying for each of the q blocks
%               on which of the M output variables it depends;
%
% OUTPUTS
% - degs_out:   1 x q cell, with degs_out{i} a 1 x Mi cell of column vectors
%               specifying the exponents of the monomials in each of the Mi
%               output variables on which block i depends, as the field ZL
%               in 'sopvar' (maximal degrees d are expanded to (0:d)');
% - degs_new:   same as degs, but with degs_new.out=0, so that 'monopmat'
%               declares operators with kernels in the input variables only;
%
% NOTES
% - Mirrors the processing of degs.out in 'monopmat' (subroutine
%   process_degs) for a single factor, as the primary monomial appears only
%   in the first factor of the tensor-PI operator (opts.otype(1)==1);

[q,M] = size(is_dep_out);
if isnumeric(degs) || iscell(degs)
    % Same degrees for input and output, as in monopmat
    degs = struct('in',{degs},'out',{degs});
elseif ~isstruct(degs) || ~isfield(degs,'in') || ~isfield(degs,'out')
    error("Options 'degs' must be specified as a struct with fields 'in' and 'out'.")
end

% Set the degrees of the output variables to 0, keeping those of the input
degs_new = degs;
degs_new.out = 0;

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
