function [prog,Pop,Qcell,alpha_list] = sopquadvar(prog,dim,vars,dom,deg,options)
% [prog,Pop,Qcell,alpha_list] = SOPQUADVAR(prog,dim,vars,dom,deg,options)
% declares a self-adjoint 'sdopvar' decision operator as a QUADRATIC FORM in
% a list of basis operators,
%
%       Pop: L_2^m[S3] -> L_2^m[S3],
%       Pop = sum_{i,j} (Z_{alpha_i})* Q_{ij} (Z_{alpha_j}),
%
% following Sec. 7 of the sopvar document. This is the constructor; it does
% not itself impose positivity. Whether the coefficient matrix Q is
% constrained is 'options.type', passed straight through to 'sosquadvar',
% exactly as 'sosquadvar' is the workhorse for 'sostools' and 'dpvar'. With
% the default type 'pos' the result satisfies Pop = Pop* >= 0, which is what
% 'possopvar' asks for; with 'sym' it is an unconstrained self-adjoint one.
%
% Each basis operator carries a single multi-index alpha and acts as
%
%       (Z_alpha x)(s) = int_theta I_alpha(s-theta) ...
%                           (I_m kron Z^alpha(theta,s)) x(theta) dtheta,
%
%       I_alpha(s) = prod_k I_{alpha_k}(s_k),
%
%       I_0(r) = delta(r),   I_1(r) = 1 for r>=0,   I_-1(r) = 1 for r<=0.
%
% Composing the adjoint pair and eliminating the integration variable theta
% by the semiseparable integral of Sec. 6.1 gives the sopvar kernel form
%
%       (Pop x)(s) = sum_gamma int_{s'} K_gamma(s,s') I_gamma(s-s') x(s') ds'
%
%       K_gamma(s,s') = (I_m kron ZL(s))' C_gamma(d) (I_m kron ZR(s'))
%
% with C_gamma(d) = unvec(A_gamma + B_gamma'*d) affine in the entries d of Q.
%
% For a single spatial variable this reproduces the L_2 -> L_2 block of
% 'poslpivar': alpha=1 is the multiplier block Z1(s), alpha=2 is the lower
% integral block Z2(eta,s), and alpha=3 is the upper integral block
% Z3(eta,s).
%
% INPUT
% - prog:   'struct' specifying an LPI/SOS program to modify;
% - dim:    scalar m, or 1x2 array [m,m], specifying the dimension of the
%           function space L_2^m. Since Pop is self-adjoint the input and
%           output dimensions must agree;
% - vars:   1 x n3 'cellstr' object specifying the names of the spatial
%           variables, or a scalar n3, in which case the variables are named
%           's1',...,'sn3'. These are the variables common to the input and
%           output space, i.e. the set S3 of the sopvar document;
% - dom:    n3 x 2 array of type 'double', with dom(k,:) = [ak,bk] the
%           domain of variable k. May be given as a single 1x2 row, in which
%           case the same domain is used for every variable;
% - deg:    degrees of the monomial bases Z^alpha(theta,s). May be given as
%             * a scalar d, in which case every basis has degree at most d
%               in each theta_k and each s_k;
%             * a 'struct' with any of the fields
%                   deg.int:    scalar or 1 x n3, maximal degree in each
%                               integration variable theta_k. Defaults to 1;
%                   deg.mult:   scalar or 1 x n3, maximal degree in each
%                               output variable s_k. Defaults to deg.int;
%                   deg.joint:  scalar, maximal total degree of theta and s
%                               combined. Defaults to no further restriction;
%                   deg.subset: 2^(2*n3) array, one maximal degree per
%                               SUBSET of the variables, over the ordering
%                               [theta_1..theta_n3, s_1..s_n3]. The bound
%                               for the subset whose members are the bits
%                               set in b sits at index 1+sum(2.^(b-1)), the
%                               convention 'poslpivar_2d/build_monoms' uses;
%                               entry 1 is unused. 'int', 'mult' and 'joint'
%                               are its singleton and full-set special
%                               cases. Optional, and needed only to
%                               reproduce a basis that per-variable plus
%                               total degrees cannot describe, as the 2D
%                               settings files do -- see
%                               'settings2possopvar';
%               For one spatial variable, (deg.int,deg.mult,deg.joint)
%               correspond to the elements of 'poslpivar's d{2};
%             * a cell array with one such scalar or struct per included
%               basis operator, in the order of the rows of alpha_list;
%           Note that Z^alpha is taken to have degree 0 in s_k whenever
%           alpha_k = 1, since the factor delta(s_k-theta_k) identifies the
%           two variables and any s_k dependence would be redundant;
% - options: (optional) 'struct' specifying other options, with fields
%   options.type      'pos' (default) declares the coefficient matrix Q
%                     positive semidefinite, so that Pop = Pop* >= 0. 'sym'
%                     declares it symmetric but otherwise unconstrained,
%                     giving a self-adjoint operator VARIABLE rather than a
%                     positive one -- the Zop of the settings files. Both
%                     values are 'sosquadvar's own vocabulary and are passed
%                     straight through. This single argument is the whole of
%                     the positivity question; the rest of this routine is
%                     assembly;
%   options.psatz     set to 1 to enforce positivity of the operator only on
%                     the domain, by including the factor
%                     g(theta) = prod_k (theta_k-ak)*(bk-theta_k) in the
%                     construction; set to 2k+1 or 2k+2 (k = 1..n3) for the
%                     linear weight of ONE face of the box,
%                     g = (theta_k-ak)/(bk-ak) or (bk-theta_k)/(bk-ak).
%                     Here k counts S3 in SORTED order -- the k-th entry of
%                     sort(vars), [ak,bk] its row of dom -- unlike dom, deg,
%                     sep and include, which follow 'vars'; for the same
%                     single space a code then names the same face as in
%                     'copquadvar' (whose registry spans all its spaces).
%                     Declare one variable per face and sum them, each at
%                     the full degree (see 'copquadvar'). Any other value
%                     is an error, including
%                     2, the ball of 'poslpivar_2d'. Defaults to 0. Rejected
%                     for type 'sym', where a nonnegative weight means
%                     nothing;
%   options.sep       logical scalar or 1 x n3 array. Where sep(k) is true,
%                     direction k is SEPARABLE: the lower and upper integral
%                     basis operators are replaced by a single full-domain
%                     integral, so the returned operator has equal lower and
%                     upper kernels in that direction. This is the R_1 = R_2
%                     form of 'poslpivar_2d's option of the same name. The
%                     admissible alpha_k are then 1 and 4, where 4 denotes
%                     the full-domain integral. Defaults to all false;
%   options.include   specification of which basis operators to include.
%                     Either an N x n3 array whose rows are the desired
%                     multi-indices alpha (entries in {1,2,3}, or {1,4} in a
%                     separable direction, with 1 the multiplier, 2 the lower
%                     integral, 3 the upper and 4 the full-domain integral),
%                     or a logical array with one entry per basis operator in
%                     their standard linear order, or a numeric vector of
%                     linear indices into that order. Defaults to all;
%
% OUTPUT
% - prog:       'struct' specifying the same program as the input, but now
%               including the decision variables defining Pop and, for type
%               'pos', a constraint enforcing Q>=0;
% - Pop:        m x m 'sdopvar' object representing a self-adjoint PI
%               operator decision variable, positive semidefinite when the
%               type is 'pos';
% - Qcell:      N x N cell array of 'cellstr' objects naming the decision
%               variables in each block Q_{ij} of the Gram matrix Q, as
%               returned by 'sosquadvar';
% - alpha_list: N x n3 array whose rows are the multi-indices alpha of the
%               included basis operators, in the order used for Qcell;
%
% NOTES
% This implements only the L_2 -> L_2 part of the construction, i.e. the
% reduction of Sec. 9.1 in which the created and lost variables S1 and S2
% are ignored. A positive operator on R^n x L_2^m[S3] additionally requires
% the R^n -> R^n, L_2 -> R^n and R^n -> L_2 blocks, which are formed from
% the same Q>=0 but need a container object to hold them.
%
% Following Sec. 9.3, the products (Z_{alpha_i})* Q_{ij} (Z_{alpha_j}) are
% formed as 'dpvar' objects via 'sosquadvar' and then converted to the
% kernel format with 'dpvar2sdvar'. This reuses existing code at the cost of
% carrying three groups of variables through the 'dpvar' arithmetic, and is
% the temporary measure the document describes; the coefficient matrices
% B_gamma could instead be built directly from vec(Q) by index arithmetic.
%
% The construction itself is 'copquadvar' on one space (10/01/2026): this
% routine is its single-space interface, in the caller's variable order.
%
% See also COPQUADVAR, POSSOPVAR, POSLPIVAR, SOSQUADVAR.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - possopvar
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
% MP, 08/22/2026: Initial coding (as 'possopvar')
% MMP, 10/01/2026: One implementation of the Sec. 9 construction: this file
%                  now checks its inputs as before, moves them to the sorted
%                  registry and calls 'copquadvar' on one space, returning
%                  the single block. The construction below the checks is
%                  deleted (see BEGIN/END). It had diverged from copquadvar
%                  five times since 09/22/2026; copquadvar held the 09/26
%                  speed-ups, this file the zero-cell skip (now in both).
%                  For sorted variables the program, Pop, Qcell and
%                  alpha_list are bit-identical to before. For unsorted ones
%                  the program equals the former one for the sorted
%                  variables, the Gram blocks list monomials in registry
%                  order, and Qcell is reordered back to the caller's
%                  monomial order so it keeps its meaning. Every entry below that
%                  describes the deleted construction (09/07, 09/11 for the
%                  cell relabelling, 09/21 moved subfunctions, 09/22 and
%                  09/30 SGCELL) no longer applies; 'process_degrees' had
%                  this file as its only caller and is deleted.
% MMP, 09/30/2026: SGCELL instrumentation removed: two counters in the gamma
%                  cell loop (MMP, 09/22/2026) and a print when the
%                  environment variable SGCELL was set, which no file in the
%                  repository sets. Diagnostic only; nothing it computed was
%                  returned.
% MMP, 09/30/2026: The reserved-suffix check and the enumeration of the
%                  basis multi-indices alpha (sopvar_implementation_notes.pdf
%                  Sec. 9.1; direction 1 fastest)
%                  call the shared 'check_reserved_names' and
%                  'multiindex_grid' (sopvar/misc/conventions) instead of
%                  inline copies of what 'copquadvar', 'lpivar_cdopvar',
%                  'eq_opts_sopvar' and 'spaces2meta_sop' also wrote. Same
%                  message, same alpha_all (test_index_maps (B) compares the
%                  verbatim inline form).
% MMP, 09/30/2026: 'lpis_sopvar/private/remap_dvars' deleted: no caller
%                  since the 09/22/2026 entry below replaced its call (kept
%                  there as a '(was)' line) by 'dvar_rows' plus a triplet
%                  scatter. The 09/21/2026 entry's list of moved subfunctions
%                  therefore no longer holds for it. Comments only.
% MMP, 09/28/2026: options.psatz is value-checked, and 2k+1 / 2k+2 weight the
%                  Gram form by one face of the box, (theta_k-a_k)/L_k or
%                  (b_k-theta_k)/L_k, k over SORTED S3. The value was taken
%                  for its truthiness, so 2 (poslpivar_2d's ball), every
%                  face code and out-of-range codes such as 9 at n3 = 2 all
%                  silently built the product weight, while 'copquadvar'
%                  (MMP, 09/27/2026) reads 3.. as faces: one options struct,
%                  two cones. Check, double conversion and face branch
%                  mirror copquadvar's; sorted k makes a code name the same
%                  face in both. Help text updated. psatz 0/1 unchanged.
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  mopquadvar -> copquadvar, posmopvar -> poscopvar.
% MMP, 09/22/2026: Defer the scatter onto the global decision variable basis
%                  until after the theta elimination. 'remap_dvars' was called
%                  before 'pack_sheets', so every block was widened from its
%                  own decision rows to the global count BEFORE the integral,
%                  and since 'pack_sheets' lays one sheet per row that put the
%                  global count on the sheet axis of every intermediate.
%                  Nothing between there and 'unpack_sheets' mixes sheets --
%                  each stays a contiguous block of columns, as
%                  'unpack_sheets' own header records -- so the widening was
%                  pure zero padding. MEASURED at two variables and degree 2,
%                  ndec = 97461: the sheet axis was 97462 columns for block
%                  pair (1,1), which needs 46, a factor of 2119, and 292386
%                  against 2190 for pair (2,3). Profiling the same case put
%                  48% of the runtime in 'int_semisep>rearrangeCoef' at 166 MB
%                  of allocation per call, with 392 GB allocated in total
%                  against a 242 MB peak -- the churn was the padding. The
%                  block now carries its own rows through the elimination and
%                  is scattered once per parameter cell, built from triplets
%                  rather than by the index-assignment into a q-row sparse
%                  that CLAUDE.md s2 forbids. Alternating A/B, two
%                  repetitions: 20.85/21.10 s becomes 5.94/6.04 s at two
%                  variables and degree 2 (3.5x), 0.18/0.19 becomes 0.15/0.11
%                  at one variable and degree 3, and the small cases are
%                  unchanged -- the saving tracks the padding factor, which
%                  grows with ndec. Coefficients are unchanged:
%                  'test_possopvar' compares the defining identity against an
%                  independent oracle and still passes 14 of 14.
% MMP, 09/21/2026: Eight subfunctions moved verbatim to 'lpis_sopvar/private'
%                  - 'process_degrees', 'process_degrees_one',
%                  'expand_caps', 'remap_dvars', 'pack_sheets',
%                  'unpack_sheets', 'lr_multiply' and 'expand_full'. Sec. 8.4
%                  asks for one implementation shared with the multi-domain
%                  'copquadvar' rather than two, and a file-local subfunction
%                  cannot be shared. A 'private' folder is visible to every
%                  function in 'lpis_sopvar' and to nothing else, so the
%                  bodies are unchanged and no name is exposed.
% MMP, 09/21/2026: Renamed from 'possopvar' and the variable type made an
%                  argument. Sec. 8.4 asks for a poscopvar over the mixed
%                  domains of a 'copvar' and observes that there should not
%                  be two implementations, just the general one with a
%                  fastlane. Extracting the constructor is the first step:
%                  the body assembles sum_ij Z_i^* Q_ij Z_j from a basis
%                  list and a coefficient matrix, and the ONLY place
%                  positivity entered was the literal 'pos' handed to
%                  'sosquadvar', so this is not a positivity routine and
%                  should not be named as one. 'possopvar' is now a wrapper
%                  fixing type to 'pos', with its behaviour unchanged, and
%                  'sym' additionally becomes available for the
%                  unconstrained operator variable the settings files call
%                  Zop. Multi-domain support, and the 'poscopvar' wrapper
%                  over it, follow separately.
% MMP, 09/21/2026: Accept 'deg.subset', one maximal degree per SUBSET of the
%                  variables, over the order the basis is built in. The
%                  previous vocabulary was one cap per variable plus a
%                  single total cap, which cannot describe the bases the 2D
%                  settings files ask for: 'poslpivar_2d/build_monoms'
%                  prunes with 2^nvars subset caps. Converting such a file
%                  could therefore only produce a covering SUPERSET,
%                  measured at 1.05x to 1.80x the basis across the six
%                  shipped 2D files, so the two builders were not the same
%                  cone at stock settings. 'int', 'mult' and 'joint' are the
%                  singleton and full-set special cases and still work
%                  exactly as before; the field is optional and empty by
%                  default, so nothing existing changes.
% MMP, 09/12/2026: Support 'options.sep', matching 'poslpivar_2d'. A
%                  separable direction replaces its lower and upper integral
%                  basis operators by one full-domain integral, encoded as
%                  alpha_k = 4, so the returned operator satisfies
%                  R_lower = R_upper there: int_a^s + int_s^b = int_a^b over
%                  a shared kernel. Implemented by expanding each 4 into its
%                  2/3 substitutions on both sides of B_i'*Q*B_j, against
%                  the same Q, and emitting each term as its own block --
%                  'plus_batch' already sums the blocks and merges their
%                  monomials, so the rest of the loop is untouched. The
%                  motivation is degree, not expressiveness: a full-domain
%                  integral has CONSTANT limits, so composing with it never
%                  substitutes a spatial variable into an integration limit,
%                  where the Volterra form inflates the monomial degree at
%                  every composition. On the 2D heat equation that inflation
%                  drives the composed derivative to per-variable degree 12
%                  against 3 for the 'opvar2d' path, and SeDuMi then fails
%                  with numerr=2 at every lambda. It also shrinks the basis
%                  from 3^n3 blocks to 2^n3 when every direction is
%                  separable, at the cost of up to 4^nsep more
%                  'int_semisep' calls per block pair.
% MMP, 09/07/2026: Collect the nblk^2 blocks and sum them once with
%                  'plus_batch' instead of accumulating pairwise. Every
%                  block is built on the same Zd, so one synchronization of
%                  decision variables and monomial bases serves all of them;
%                  accumulating pairwise re-synchronized the growing
%                  accumulator on every addition, which profiled at 32% of
%                  this routine at three spatial variables (729 additions).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % Matrix dimension
if isscalar(dim)
    m = dim;
elseif numel(dim)==2
    if dim(1)~=dim(2)
        error("A self-adjoint operator must have equal input and output dimensions.")
    end
    m = dim(1);
else
    error("Matrix dimension of the operator should be specified as a scalar or 1x2 array.")
end
if m<1 || m~=round(m)
    error("Matrix dimension of the operator should be a positive integer.")
end

% % % Spatial variables
if isnumeric(vars) && isscalar(vars)
    n3 = vars;
    vars = cell(1,n3);
    for k=1:n3
        vars{k} = ['s',num2str(k)];
    end
elseif isa(vars,'polynomial')
    vars = vars.varname(:)';
elseif ischar(vars)
    vars = {vars};
end
if ~iscellstr(vars)
    error("Spatial variables should be specified as a 1 x n3 'cellstr' object.")
end
vars = vars(:)';        n3 = numel(vars);
if numel(unique(vars))~=n3
    error("Spatial variable names should be unique.")
end
% The parameter cell of an sopvar/sdopvar is indexed over S3 in SORTED      % MMP, 09/11/2026
% order: the class recovers a direction from a cell subscript with          % MMP, 09/11/2026
% intersect(vars.in,vars.out), which sorts (see canonicalize_multiplier,    % MMP, 09/11/2026
% apply_sopvar and lpi_eq_sdopvar). Everything below is built over 'vars'   % MMP, 09/11/2026
% in the CALLER's order, because 'dom' and 'alpha_list' have to stay        % MMP, 09/11/2026
% consistent with each other for 'int_semisep'. This records the            % MMP, 09/11/2026
% permutation so the finished cell can be relabelled onto sorted order      % MMP, 09/11/2026
% just before the constructor is called; sorted variable j is              % MMP, 09/11/2026
% vars{ord_S3(j)}.                                                          % MMP, 09/11/2026
[~,ord_S3] = sort(vars);                                                    % MMP, 09/11/2026
% (10/01/2026: the comment above describes the deleted construction; ord_S3 % MMP, 10/01/2026
% now maps the caller's order onto copquadvar's sorted registry, below.)    % MMP, 10/01/2026

% Declare names for the integration variable theta and the input dummy
% variable s'. Neither appears in the returned object: the dummy variables
% of an sdopvar are implicit in vars.in.
% These suffixes are reserved outright, rather than only where they would
% actually collide, so that the set of admissible names does not depend on
% which other variables happen to be present.
% vars_int = strcat(vars,'_int');                                           % MMP, 10/01/2026 (was)
% vars_dum = strcat(vars,'_dum');                                           % MMP, 10/01/2026 (was)
% copquadvar builds both name lists now; the reservation stays here.        % MMP, 10/01/2026
% is_reserved = ~cellfun(@isempty,regexp(vars,'_(int|dum)$','once'));       % MMP, 09/30/2026 (was)
% if any(is_reserved)                                                       % MMP, 09/30/2026 (was)
%     error("Spatial variable names may not end in '_int' or '_dum'; "...
%           +"'"+string(vars{find(is_reserved,1)})+"' does.")               % MMP, 09/30/2026 (was)
% end                                                                       % MMP, 09/30/2026 (was)
check_reserved_names(vars);     % shared check, same message                % MMP, 09/30/2026

% % % Domain
if isempty(dom) && n3==0
    dom = zeros(0,2);
end
if size(dom,2)~=2
    error("Domains should be specified as an n3 x 2 array.")
end
if size(dom,1)==1 && n3~=1
    dom = repmat(dom,n3,1);
elseif size(dom,1)~=n3
    error("Domains should be specified as an n3 x 2 array for n3 spatial variables.")
end
if any(dom(:,2)<=dom(:,1))
    error("Each domain should satisfy dom(k,1) < dom(k,2).")
end

% % % Options
if nargin<6 || isempty(options)
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if ~isfield(options,'psatz') || isempty(options.psatz)
    options.psatz = 0;
end
% 0, 1 and the face codes 2k+1 / 2k+2, k over SORTED S3 as in 'copquadvar'; % MMP, 09/28/2026
% checked, not taken for truthiness, which let 2 (poslpivar_2d's ball) and  % MMP, 09/28/2026
% every face code build the product while 'copquadvar' reads 3.. as faces.  % MMP, 09/28/2026
if ~isscalar(options.psatz) || ~ismember(options.psatz,[0,1,3:2*n3+2])      % MMP, 09/28/2026
    error("'psatz' should be 0, 1, or a face code 2k+1 / 2k+2 with k "...
          +"in 1..n3, here n3 = "+num2str(n3)+"; 2 is "...
          +"poslpivar_2d's ball, not the box.")                             % MMP, 09/28/2026
end                                                                         % MMP, 09/28/2026
% Double from here on: an integer class passes the check above, and then    % MMP, 09/28/2026
% floor((psatz-1)/2) rounds in integer arithmetic and picks the wrong face. % MMP, 09/28/2026
options.psatz = double(options.psatz);                                      % MMP, 09/28/2026
% % % Variable type. This is the ONLY place positivity enters: everything    % MMP, 09/21/2026
% else here assembles sum_ij Z_i^* Q_ij Z_j from the basis list and the      % MMP, 09/21/2026
% coefficient matrix, and whether Q is constrained is one argument to        % MMP, 09/21/2026
% 'sosquadvar'. 'pos' reproduces what 'possopvar' has always declared;       % MMP, 09/21/2026
% 'sym' gives an unconstrained self-adjoint operator variable, the Zop of   % MMP, 09/21/2026
% the settings files. Both names are 'sosquadvar' own vocabulary.           % MMP, 09/21/2026
vartype = 'pos';                                                            % MMP, 09/21/2026
if isfield(options,'type') && ~isempty(options.type)                        % MMP, 09/21/2026
    vartype = char(options.type);                                           % MMP, 09/21/2026
    if ~ismember(vartype,{'pos','sym'})                                    % MMP, 09/21/2026
        error("'type' should be 'pos' or 'sym'.")                          % MMP, 09/21/2026
    end                                                                     % MMP, 09/21/2026
end                                                                         % MMP, 09/21/2026
if strcmp(vartype,'sym') && options.psatz                                  % MMP, 09/21/2026
    % A psatz weight multiplies the form by g(theta) >= 0 on the domain,     % MMP, 09/21/2026
    % which only means anything when Q is constrained.                       % MMP, 09/21/2026
    error("'psatz' has no meaning for a 'sym' variable; use type 'pos'.")  % MMP, 09/21/2026
end                                                                         % MMP, 09/21/2026
% % % Separable directions                                                  % MMP, 09/12/2026
% 'sep(k)' replaces the lower and upper integral basis operators in         % MMP, 09/12/2026
% direction k by a single FULL-DOMAIN integral, encoded as alpha_k = 4.     % MMP, 09/12/2026
% Since int_a^s + int_s^b = int_a^b over one shared kernel, the operator    % MMP, 09/12/2026
% then has equal lower and upper kernels in that direction: the R_1 = R_2   % MMP, 09/12/2026
% form 'poslpivar_2d' calls 'sep'. Its limits are constant, so composing    % MMP, 09/12/2026
% with it never substitutes a spatial variable into an integration limit    % MMP, 09/12/2026
% and never inflates the monomial degree, which the Volterra form does at   % MMP, 09/12/2026
% every composition.                                                        % MMP, 09/12/2026
if ~isfield(options,'sep') || isempty(options.sep)                          % MMP, 09/12/2026
    sep = false(1,n3);                                                      % MMP, 09/12/2026
else                                                                        % MMP, 09/12/2026
    sep = logical(reshape(options.sep,1,[]));                               % MMP, 09/12/2026
    if isscalar(sep)                                                        % MMP, 09/12/2026
        sep = repmat(sep,1,n3);                                             % MMP, 09/12/2026
    elseif numel(sep)~=n3                                                   % MMP, 09/12/2026
        error("'sep' should be a scalar or have one entry per variable.")   % MMP, 09/12/2026
    end                                                                     % MMP, 09/12/2026
end                                                                         % MMP, 09/12/2026

% % % Basis operators to include
%   1 <-> multiplier (delta),  2 <-> lower integral,  3 <-> upper integral,
%   4 <-> full-domain integral, used exactly where 'sep' is set
% Direction 1 varies fastest, the layout 'fliplr(dec2base(...))' produced.  % MMP, 09/12/2026
vals = cell(1,n3);                                                          % MMP, 09/12/2026
for k = 1:n3                                                                % MMP, 09/12/2026
    if sep(k), vals{k} = [1,4]; else, vals{k} = [1,2,3]; end                % MMP, 09/12/2026
end                                                                         % MMP, 09/12/2026
% szv = cellfun(@numel,vals);                                               % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
% nall = prod([szv,1]);                                                     % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
% alpha_all = zeros(nall,n3);                                               % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
% rep = 1;                                                                  % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
% for k = 1:n3                                                              % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
%     col = reshape(repmat(vals{k},rep,1),[],1);                            % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
%     alpha_all(:,k) = repmat(col,nall/(rep*szv(k)),1);                     % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
%     rep = rep*szv(k);                                                     % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
% end                                                                       % MMP, 09/12/2026 % MMP, 09/30/2026 (was)
alpha_all = multiindex_grid(vals);  % the same loop, shared (conventions/)  % MMP, 09/30/2026
if ~isfield(options,'include') || isempty(options.include)
    alpha_list = alpha_all;
else
    incl = options.include;
    if islogical(incl)
        if numel(incl)~=size(alpha_all,1)
            error("A logical 'include' has one entry per basis operator.")  % MMP, 09/12/2026
        end
        alpha_list = alpha_all(incl(:),:);
    elseif n3>0 && size(incl,2)==n3
        % A separable direction admits only 1 or 4 and a non-separable one  % MMP, 09/12/2026
        % only 1, 2 or 3, so validate against the generated list rather     % MMP, 09/12/2026
        % than a fixed numeric range.                                       % MMP, 09/12/2026
        if ~all(ismember(incl,alpha_all,'rows'))                            % MMP, 09/12/2026
            error("'include' contains an inadmissible multi-index.")        % MMP, 09/12/2026
        end
        alpha_list = incl;
    else
        if any(incl(:)<1) || any(incl(:)>size(alpha_all,1)) ...
                || any(incl(:)~=round(incl(:)))
            error("Linear indices in 'include' are out of range.")          % MMP, 09/12/2026
        end
        alpha_list = alpha_all(incl(:),:);
    end
end
% Note that for n3=0 there is exactly one basis operator, indexed by the
% empty multi-index, so the row count rather than 'isempty' must be tested.
if size(alpha_list,1)==0
    error("At least one basis operator must be included.")
end
if size(unique(alpha_list,'rows'),1)~=size(alpha_list,1)
    error("Multi-indices in 'include' should be distinct.")
end
nblk = size(alpha_list,1);

% BEGIN MMP, 10/01/2026: the construction is 'copquadvar' on one space (spec
% sec. 9; the Overleaf sec. 8.4 asks for one implementation, the single
% space being a special case). Deleted here, old lines 497-748: the degree
% processing ('process_degrees'), the monomial bases, the 'sosquadvar'
% call, the psatz weight and the block-pair loop with its 'int_semisep',
% 'unpack_sheets', 'lr_multiply' and 'plus_batch' assembly, stamped MP
% 08/22/2026 and MMP 09/07, 09/11, 09/21, 09/22, 09/25, 09/28, 09/30/2026;
% git history holds them. The checks above are unchanged, so the accepted
% inputs and the error messages are as before.
%
% copquadvar works over the SORTED registry; registry variable j is
% vars{ord_S3(j)}. 'dom' goes in paired with the names, so it needs no
% permutation; 'sep', the degree caps and the multi-indices are permuted.
% The basis list built above goes in as 'include', which copquadvar keeps
% in the given order, so alpha_list, Qcell and a per-basis 'deg' cell keep
% this routine's order.
o = options;                                                                % MMP, 10/01/2026
o.sep = sep(ord_S3);                                                        % MMP, 10/01/2026
o.include = alpha_list(:,ord_S3);                                           % MMP, 10/01/2026
if iscell(deg)                                                              % MMP, 10/01/2026
    if numel(deg)~=nblk                                                     % MMP, 10/01/2026
        error("A cell 'deg' should have one entry per included basis operator.") % MMP, 10/01/2026
    end                                                                     % MMP, 10/01/2026
    degc = {cellfun(@(d) registry_deg(d,ord_S3,n3),reshape(deg,1,[]), ...
                    'UniformOutput',false)};    % one space, one entry per basis % MMP, 10/01/2026
else                                                                        % MMP, 10/01/2026
    degc = registry_deg(deg,ord_S3,n3);                                     % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
[prog,Pc,Qcell] = copquadvar(prog,m,{vars},struct('vars',{vars},'dom',dom),degc,o); % MMP, 10/01/2026
% Unsorted variables: copquadvar lists each basis operator's monomials over % MMP, 10/01/2026
% the registry, this routine listed them over 'vars'. Rebuild both exponent % MMP, 10/01/2026
% grids with the helpers copquadvar uses and reorder the rows and columns   % MMP, 10/01/2026
% of each Qcell block (sosquadvar: component outer, monomial inner), so that % MMP, 10/01/2026
% Qcell keeps its meaning. The program and Pop are not affected.            % MMP, 10/01/2026
if ~isequal(reshape(ord_S3,1,[]),1:n3)                                      % MMP, 10/01/2026
    Rw = cell(1,nblk);                                                      % MMP, 10/01/2026
    for i = 1:nblk                                                          % MMP, 10/01/2026
        spec = deg;     if iscell(deg),  spec = deg{i};  end                % MMP, 10/01/2026
        dc = process_degrees_one(spec,n3);                                  % MMP, 10/01/2026
        dr = process_degrees_one(registry_deg(spec,ord_S3,n3),n3);          % MMP, 10/01/2026
        mc = dc.mult;   mc(alpha_list(i,:)==1) = 0;                         % MMP, 10/01/2026
        mr = dr.mult;   mr(alpha_list(i,ord_S3)==1) = 0;                    % MMP, 10/01/2026
        Ec = build_exponent_grid([dc.int,mc],dc.joint,dc.subset);           % MMP, 10/01/2026
        Er = build_exponent_grid([dr.int,mr],dr.joint,dr.subset);           % MMP, 10/01/2026
        [tf,p] = ismember(Ec(:,[ord_S3,n3+ord_S3]),Er,'rows');              % MMP, 10/01/2026
        if ~all(tf) || numel(p)~=size(Er,1)                                 % MMP, 10/01/2026
            error("Internal error: the monomial bases do not correspond.")  % MMP, 10/01/2026
        end                                                                 % MMP, 10/01/2026
        Rw{i} = reshape(p(:)+size(Er,1)*(0:m-1),[],1);                      % MMP, 10/01/2026
    end                                                                     % MMP, 10/01/2026
    for i = 1:nblk                                                          % MMP, 10/01/2026
        for j = 1:nblk                                                      % MMP, 10/01/2026
            Qcell{i,j} = Qcell{i,j}(Rw{i},Rw{j});                           % MMP, 10/01/2026
        end                                                                 % MMP, 10/01/2026
    end                                                                     % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
Pop = Pc.C{1,1};                                                           % MMP, 10/01/2026
% END MMP, 10/01/2026

end


function d = registry_deg(d,ord,n3)
% One degree specification moved from the caller's variable order to the
% sorted registry: the per-variable caps by ORD, and the 'subset' array by
% relabelling the bits of its index over [theta_1..theta_n3, s_1..s_n3]
% (see the help: the subset with bits b sits at 1+sum(2.^(b-1))). A scalar
% or a scalar cap is the same in every order.
%
% Initial coding MMP, 10/01/2026
if ~isstruct(d) || n3<2                                                     % MMP, 10/01/2026
    return                                                                  % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
for f = {'int','mult'}                                                      % MMP, 10/01/2026
    if isfield(d,f{1}) && numel(d.(f{1}))==n3                               % MMP, 10/01/2026
        d.(f{1}) = d.(f{1})(ord);                                           % MMP, 10/01/2026
    end                                                                     % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
if isfield(d,'subset') && numel(d.subset)==2^(2*n3)                         % MMP, 10/01/2026
    r(ord) = 1:n3;                              % registry position of vars{k} % MMP, 10/01/2026
    nb = [r, n3+r];                             % caller bit b -> registry bit % MMP, 10/01/2026
    bits = fliplr(dec2bin(0:2^(2*n3)-1,2*n3)=='1');     % column b = bit b  % MMP, 10/01/2026
    to = bits*(2.^(nb(:)-1)) + 1;               % registry index of each entry % MMP, 10/01/2026
    s = d.subset;                                                           % MMP, 10/01/2026
    s(to) = d.subset(:);                                                    % MMP, 10/01/2026
    d.subset = s;                                                           % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
end                                                                         % MMP, 10/01/2026
