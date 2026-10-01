function [prog,Pop,Qcell,basis_list] = copquadvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,Qcell,basis_list] = COPQUADVAR(prog,dims,spaces,dom,deg,options)
% declares a self-adjoint 'cdopvar' decision operator on a CONCATENATION OF
% MIXED SPACES as a quadratic form in a list of basis operators, following
% Sec. 8.4 of the sopvar document,
%
%       Pop: X -> X,     X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M],
%
%       Pop = sum_{k,l} sum_{i,j} (Z_{alpha_i^k})* Q_{(ki),(lj)} (Z_{alpha_j^l})
%
% so that block (k,l) of the container is the (k,l) group of that sum. There
% is ONE coefficient matrix Q spanning every (space, multi-index) pair, which
% is what makes the whole container positive rather than only its diagonal.
% An empty space is R^{m_k}, so this covers the mixed finite/infinite case
% that every PIE has and that 'sopquadvar' explicitly does not: see its NOTES.
%
% Each basis operator carries a multi-index over its OWN space's variables
% and acts into a common auxiliary space L_2[theta] over all the variables,
%
%       (Z_alpha x)(theta) = int_{s in dom(s^k)} I_alpha(theta-s) ...
%                               (I_{m_k} kron Z^alpha(theta,s)) x(s) ds,
%
%       I_alpha(r) = prod_d I_{alpha_d}(r_d),
%
%       I_0(r) = delta(r),  I_1(r) = 1 for r>=0,  I_-1(r) = 1 for r<=0,
%
% with the product indicator taken over the variables of s^k only: a variable
% the space does not have gets no indicator, and the basis depends on its
% theta polynomially. Pairing two such operators integrates theta out, and
% each direction d falls into one of four cases,
%
%   d in s^k and s^l   two-sided, the semiseparable integral of Sec. 6.1;
%                      d is a shared variable S3 of block (k,l)
%   d in s^k only      one-sided, int_dom I_beta(s_d-theta_d) . dtheta_d;
%                      d is a created variable S2, a free multiplier
%   d in s^l only      one-sided, int_dom I_alpha(theta_d-s'_d) . dtheta_d;
%                      d is a lost variable S1, integrated over its domain
%   d in neither       int_dom . dtheta_d, a plain definite integral
%
% Only the first needs 'int_semisep'; the other three are antiderivatives of
% a monomial, evaluated at a limit that is either a spatial variable or an
% endpoint. They are applied first, so 'int_semisep' sees only the shared
% directions of the pair.
%
% For a single space this reproduces 'sopquadvar' exactly - every direction
% is then shared - and for the two spaces (R^n, L_2^m[s]) it reproduces the
% four blocks of 'poslpivar'.
%
% INPUT
% - prog:   'struct' specifying an LPI/SOS program to modify;
% - dims:   M x 1 array of positive integers, dims(k) the dimension m_k of
%           space k. A scalar is used for every space. Since Pop is
%           self-adjoint the input and output spaces coincide;
% - spaces: 1 x M cell of 'cellstr' objects, spaces{k} naming the spatial
%           variables of space k. An EMPTY entry is the finite-dimensional
%           space R^{m_k}. A plain 'cellstr' is read as a single space, so
%           that copquadvar(prog,m,{'s1'},...) means what it looks like;
%           several spaces must therefore be nested, as in {{},{'s1'}}. May
%           also be given as a 'struct' with field 'out';
% - dom:    nv x 2 array, dom(d,:) = [ad,bd] the domain of the d-th variable
%           of the registry, which is the SORTED union of all the spaces. May
%           be given as a single 1x2 row, used for every variable, or as a
%           'struct' with fields 'vars' and 'dom' pairing each name with its
%           interval, which removes the need to know the registry's order;
% - deg:    degrees of the monomial bases Z^alpha(theta,s), in the vocabulary
%           of 'sopquadvar': a scalar, or a 'struct' with fields 'int',
%           'mult', 'joint' and 'subset'. 'int' and 'mult' are given over the
%           whole registry, and space k uses the 'mult' entries of its own
%           variables; the basis of space k is therefore built over
%           [theta_1..theta_nv, s^k], whose length varies with the space, so
%           'subset' - whose indexing depends on that length - is accepted
%           only for a space holding every variable. May also be a cell with
%           one such specification per SPACE, each of which may itself be a
%           cell with one per basis operator of that space;
%           the fields cap degrees of the factor Z^alpha(theta,s), not of   % MMP, 09/30/2026
%           the kernel of Pop:                                              % MMP, 09/30/2026
%             'int'    scalar or 1 x nv, the degree in each integration     % MMP, 09/30/2026
%                      variable theta_d. Defaults to 1;                     % MMP, 09/30/2026
%             'mult'   scalar or 1 x nv, the degree in each output          % MMP, 09/30/2026
%                      variable s_d, set to 0 where alpha_d is a            % MMP, 09/30/2026
%                      multiplier. Defaults to 'int';                       % MMP, 09/30/2026
%             'joint'  scalar, the total degree in (theta,s^k). Defaults    % MMP, 09/30/2026
%                      to no further cap;                                   % MMP, 09/30/2026
%             'subset' 2^(2*nv) entries, one cap per subset of              % MMP, 09/30/2026
%                      [theta_1..theta_nv, s^k], as in 'sopquadvar';        % MMP, 09/30/2026
%           'lpivar_cdopvar' names its fields 'int' and 'mult' too, with    % MMP, 09/30/2026
%           other meanings: kernel degrees per variable role, 'mult' the    % MMP, 09/30/2026
%           left degree in a multiplier direction and 'int' the [left       % MMP, 09/30/2026
%           right] degrees in an integral direction. A 'deg' struct         % MMP, 09/30/2026
%           therefore does not carry over between the two routines;         % MMP, 09/30/2026
% - options: (optional) 'struct' specifying other options, with fields
%   options.type      'pos' (default) declares Q positive semidefinite, so
%                     that Pop = Pop* >= 0; 'sym' declares it symmetric but
%                     otherwise unconstrained. Passed to 'sosquadvar', whose
%                     vocabulary this is, and the only place positivity
%                     enters;
%   options.psatz     set to 1 to enforce positivity only on the domain, by
%                     including the factor g(theta) = prod_d
%                     (theta_d-ad)*(bd-theta_d). Defaults to 0, and rejected
%                     for type 'sym', where a nonnegative weight means
%                     nothing;
%                     set to 2d+1 or 2d+2 (d = 1..nv, over the registry) for
%                     the linear weight of ONE face of the box,
%                     g = (theta_d-ad)/(bd-ad) or (bd-theta_d)/(bd-ad). Declare
%                     one variable per face and sum them, each at the full
%                     degree of the base variable (measured in
%                     'poslpivar_2d', CC 09/23/2026: face terms at eq_deg-1
%                     destroyed the certificate). In 2-D with sorted names
%                     these are 'poslpivar_2d' psatz 3-6;
%   options.sep       logical scalar or 1 x nv array over the registry. Where
%                     sep(d) is true the lower and upper integral basis
%                     operators in direction d are replaced by one
%                     full-domain integral, encoded as alpha_d = 4, exactly
%                     as in 'sopquadvar';
%   options.include   1 x M cell, entry k selecting the basis operators of
%                     space k in any of the formats 'sopquadvar' accepts: an
%                     array of multi-indices over space k's own variables, a
%                     logical mask over its basis operators, or linear
%                     indices into them. For a single space the cell may be
%                     omitted. Defaults to all;
%
% OUTPUT
% - prog:       the same program, now carrying the decision variables of Pop
%               and, for type 'pos', the constraint Q>=0;
% - Pop:        M x M 'cdopvar' container representing a self-adjoint PI
%               operator decision variable on X, positive semidefinite when
%               the type is 'pos';
% - Qcell:      N x N cell of 'cellstr' objects naming the decision variables
%               of each block of Q, as returned by 'sosquadvar', where N is
%               the total number of basis operators over all spaces;
% - basis_list: N x (1+nv) array, row c describing basis operator c: column 1
%               is its space index k, and column 1+d its multi-index entry
%               in direction d, or 0 where variable d is not in s^k. The rows
%               are ordered by space and match Qcell;
%
% NOTES
% The container's blocks are assembled from one Q, so the object returned is
% positive as a whole; constraining only the diagonal blocks would be a
% strictly weaker condition and is not what this builds.
%
% Sec. 8.4 observes that there should not be two implementations of this
% construction, just the general one with a fastlane for the single-space
% case. 'sopquadvar' is that fastlane and still has its own pair loop; the
% input processing, degree handling and coefficient bookkeeping are already
% shared through 'lpis_sopvar/private'.
% The space, domain and registry parse is 'parse_copvar_spaces', shared     % MMP, 09/30/2026
% with 'lpivar_cdopvar'.                                                    % MMP, 09/30/2026
%
% See also SOPQUADVAR, POSCOPVAR, POSSOPVAR, POSLPIVAR, SOSQUADVAR,
% INT_SEMISEP, CDOPVAR, LPI_EQ_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - copquadvar
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
% Initial coding MMP, 09/21/2026
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Renamed here:
%                  mopquadvar -> copquadvar, test_posmopvar -> test_poscopvar.
%                  File was 'mopquadvar.m'.
% MMP, 09/22/2026: Defer the scatter onto the global decision variable basis
%                  until after the elimination, the same change 'sopquadvar'
%                  took on this date and for the same reason: 'remap_dvars'
%                  before 'pack_sheets' put the global decision count on the
%                  sheet axis of every intermediate, where each block needs
%                  only its own. 'int_onesided' now receives nloc rather than
%                  ndec as its sheet count, so its triplet remapping works on
%                  the narrow axis too. Alternating A/B, two repetitions, on
%                  the four spaces R x L2[x] x L2[y] x L2[x,y]: at degree 2
%                  and ndec = 166176, 83.2/78.1 s becomes 14.3/13.5 s (5.8x);
%                  at degree 1 and ndec = 10440, 5.16/4.87 becomes 4.23/4.03.
%                  As in 'sopquadvar' the saving tracks the padding factor and
%                  so grows with ndec. Coefficients unchanged: 'test_poscopvar'
%                  checks the defining identity against an independent oracle
%                  and still passes 19 of 19.
% MMP, 09/21/2026: Force every 'repelem' result in 'int_onesided' to a column.
%                   'repelem' of a scalar returns a ROW, so with exactly one
%                   nonzero in the block's coefficient matrix the index
%                   arithmetic implicitly expanded to a matrix, each entry was
%                   emitted 'cc' times at the same position, and 'sparse'
%                   summed the duplicates - a block scaled by the one-sided
%                   fan-out, silently. Only reachable at monomial degree 0,
%                   which is why every earlier test missed it, and invisible
%                   to a span comparison since scaling a generator does not
%                   change the span.
% MMP, 09/21/2026: Validate the VALUE of options.psatz instead of its
%                   truthiness. psatz=2 was accepted and silently built the
%                   box, where 'poslpivar_2d' builds a ball - a different
%                   cone with no warning.
% MMP, 09/21/2026: Bound the per-pair term preallocation by the separable
%                   SHARED directions rather than by all shared directions.
%                   With the default sep it over-allocated an empty cell
%                   array by 4^n3 per block pair.
% MMP, 09/26/2026: Form kron(Rmat.',Lmat).' once per (iL,iR) pair and pass
%                  it to 'lr_multiply', which formed it once per gamma cell
%                  with the same Lmat and Rmat (3^n3 cells). Measured on a
%                  2-D stability build: kron was 3.2 of 11.7 s. Same result;
%                  cost unchanged in q, one kron per pair instead of 3^n3.
% MMP, 09/26/2026: Call site only: the per-gamma unpack_sheets, lr_multiply,
%                  find and scatter onto Zd are now one 'unpack_sheets'
%                  call, which applies Lmat and Rmat as an index remap of
%                  the triplets of Cgam{q}. SUPERSEDES the entry above: no
%                  Kronecker product is formed any more. With the
%                  vectorized 'dvar_rows' of the same date, alternating A/B
%                  medians of whole cx_exec builds: 2-D Hinf (io2, light,
%                  ndec 7.6e5) 55.0 -> 40.0 s; 2-D stability 9.4 -> 6.7 s
%                  light, 67.3 -> 45.3 s heavy (ndec 9.1e5); nv=3 poscopvar
%                  3.2 -> 2.5 s. Programs bit-identical in every field.
%                  What remains of this loop is 'find' on Cgam, whose
%                  column count 'int_semisep' sets at (1+nloc)*g2b*NR3.
% MMP, 09/26/2026: Call int_semisep with layout 'triplets', so Cgam{q}
%                  reaches 'unpack_sheets' as (i,j,v) and the
%                  (1+nloc)*g2b*NR3-column matrix above is never built nor
%                  scanned (4.6 + 3.6 s of a 19.7 s 2-D Hinf build). And Zd
%                  is concatenated once instead of grown per Qcell block.
%                  Programs bit-identical; measured numbers in int_semisep.
% MMP, 09/26/2026: plus_batch(...,'shared_Zd'): every term of a block is
%                  built on this call's one Zd, so its per-operand O(q) list
%                  comparison is redundant (plus_batch 1.13 -> 0.33 s, 2-D
%                  Hinf). Pblk.dvarname passed without a cellstr(string())
%                  round trip when already a cellstr. Bit-identical.
% MMP, 09/27/2026: options.psatz = 2d+1 / 2d+2 weights the Gram form by the
%                  normalised linear generator (theta_d-a_d)/L_d or
%                  (b_d-theta_d)/L_d of one face of the box, d indexing the
%                  registry: the N-D extension of 'poslpivar_2d' psatz 3-6
%                  (identical codes in 2-D with sorted names). The product
%                  weight (psatz=1) vanishes on every face, so it cannot
%                  supply what a certificate needs near one face; measured on
%                  the heatNd benchmark (u_t = Lap u + r u), 2-D: product
%                  certifies nothing, the 4 faces reach 0.9875 kappa*.
%                  Nonnegative on the box by construction, so no weight can
%                  make the variable indefinite. psatz 0/1 unchanged.
% MMP, 09/29/2026: Skip an all-zero gamma cell before 'unpack_sheets' and
%                  write its zeros directly: the skip 'sopquadvar' took on
%                  09/22/2026, which this copy lacked. At most 2^n3 of the
%                  3^n3 cells can be nonzero (Sec. 6.1 key table), so
%                  79884 of 87696 unpack calls in the 3-D heatNd build
%                  (w3) were on an empty cell. Measured unprofiled against
%                  the HEAD copy on the same arguments: 52 -> 15 us per
%                  empty cell (median, 60 captured w3 sets), ~3.0 s per w3
%                  build; the 14 w3 calls 43.5 -> 39.1 s. Cost in q
%                  unchanged: the zeros written are those the unpack
%                  returned. Programs bit-identical (w1-w3, 91 SDP
%                  leaves). An empty B stored as [] was rejected:
%                  'plus_batch' then errors ('Bsum + []').
% MMP, 09/30/2026: Shared helpers (sopvar/misc/conventions) in place of
%                  two local copies, each of which had hand-written twins
%                  elsewhere that could drift apart silently: (1) spaces,
%                  registry, dims, reserved names and domain are read by
%                  'parse_copvar_spaces', the parser shared with
%                  'lpivar_cdopvar'; the self-adjoint check stays here,
%                  after the parse. (2) 'multiindex_grid' replaces the
%                  local 'alpha_grid' (same body, direction 1 fastest,
%                  the order 'sopquadvar' and 'eq_opts_sopvar' use).
%                  'int_onesided' keeps its inline copy of 'kron_split',
%                  which the helper calls would slow (see its comment).
%                  Same values on every input the old parse accepted:
%                  'test_space_parser' compares the two on 120 random
%                  inputs, and the programs are bit-identical (w1-w3,
%                  91 SDP leaves). On other inputs the parser's NOTES
%                  list the differences: six messages reworded; the
%                  self-adjoint error now follows the other input
%                  errors; a struct ARRAY 'spaces' is refused (its first
%                  element was used); '' is one R^q space (a variable
%                  named '' was declared); struct dims and [] spaces are
%                  accepted. Time unchanged: the 14 w3 calls, run on the
%                  same arguments within one build, unprofiled, take
%                  34.99 s before and 34.88 s after (a second copy of the
%                  old file: 34.57 s).
% MMP, 09/30/2026: Split the primary function into local functions, one
%                  per phase of the construction, bodies moved unchanged:
%                  parse_inputs, parse_options, basis_operators,
%                  degree_list, monomial_bases, common_dvars, psatz_weight
%                  and block_terms, with, inside block_terms, gamma_params
%                  (the parameter cell of one term over the gamma cells of
%                  Sec. 6.1) and by_space (four identical scatter loops).
%                  Two moved lines change, marked (was): the nargin test
%                  on options moves to the primary, and the unused nlv
%                  goes. Why: the basis operators Z_alpha (Sec. 9.1) and the
%                  pairs Z_i* Q_ij Z_j (Sec. 9.2) sat behind about 260
%                  lines of input normalization, and the computation of one
%                  term 6 loops deep in one function of McCabe 101.
%                  Primary: 354 -> 24 code lines, McCabe 101 -> 5, nesting
%                  depth 8 -> 2; largest local McCabe 20 (parse_options),
%                  deepest nesting 4 (block_terms). The q-length Zd and
%                  the Nprod cell reach block_terms by reference and are
%                  only read. Programs bit-identical (w1-w3, 91 SDP
%                  leaves); the 14 w3 calls on the same arguments take
%                  34.94 s against 34.99 s before (unprofiled, one build).
%                  Also documented, in the help and in degree_list, what
%                  each 'deg' field means here, and that 'lpivar_cdopvar'
%                  gives 'int' and 'mult' other meanings; no option changed.


% % % BEGIN change MMP, 09/30/2026 (split): the primary runs the            % MMP, 09/30/2026
% construction of Sec. 9 phase by phase, one call each: the inputs; the     % MMP, 09/30/2026
% basis operators Z_alpha (Sec. 9.1) and their degree caps; their           % MMP, 09/30/2026
% monomial bases; one Gram Q through 'sosquadvar' (Sec. 9.3); the pairs     % MMP, 09/30/2026
% Z_i* Q_ij Z_j (Sec. 9.2), block by block. Each phase was a section of     % MMP, 09/30/2026
% this function and is now the local function below named after it, its     % MMP, 09/30/2026
% lines moved unchanged with their markers. New, and marked: the lines      % MMP, 09/30/2026
% passing variables in and out, and in 'block_terms' the calls to           % MMP, 09/30/2026
% 'gamma_params' and 'by_space'.                                            % MMP, 09/30/2026

if nargin<5
    error("Not enough input arguments.")
end
if nargin<6,    options = [];   end                                         % MMP, 09/30/2026

S = parse_inputs(dims,spaces,dom);                                          % MMP, 09/30/2026
[options,vartype,S.sep,incl] = parse_options(options,S.nv,S.M);             % MMP, 09/30/2026
[alpha_sp,nb,base,basis_list,bsp,bix] = basis_operators(S,incl);            % MMP, 09/30/2026
deg_list = degree_list(deg,S,nb,bsp,bix);                                   % MMP, 09/30/2026
[Z1c,Z2c,mdim] = monomial_bases(S,alpha_sp,deg_list,bsp,bix);               % MMP, 09/30/2026

% One call, one Gram, spanning every (space, multi-index) pair: this is what
% couples the blocks of the container and makes Pop positive as a whole.
[prog,Nprod,Qcell] = sosquadvar(prog,Z1c,Z2c,mdim,mdim,vartype);

% What the block loop reads. No callee assigns into Zd, Nprod or prog, so   % MMP, 09/30/2026
% none of the q-length data is copied (copy-on-write).                      % MMP, 09/30/2026
B = struct('alpha_sp',{alpha_sp},'nb',nb,'base',base,'Nprod',{Nprod});      % MMP, 09/30/2026
[B.Zd,B.dmap,B.ndec] = common_dvars(Qcell);                                 % MMP, 09/30/2026
B.gfun = psatz_weight(options,S);                                           % MMP, 09/30/2026

M = S.M;                                                                    % MMP, 09/30/2026
Cblk = cell(M,M);
for k = 1:M
  for l = 1:M
    [terms,nt] = block_terms(k,l,S,B);                                      % MMP, 09/30/2026
%   Cblk{k,l} = plus_batch(terms{1:nt});                                    % MMP, 09/26/2026 (was)
    % Every term was built on this call's one Zd (the sdopvar constructor   % MMP, 09/26/2026
    % stores it unchanged), so plus_batch may skip its O(q) per-operand     % MMP, 09/26/2026
    % list comparison: 136 comparisons of 3.8e5 names in 2-D Hinf.          % MMP, 09/26/2026
    Cblk{k,l} = plus_batch(terms{1:nt},'shared_Zd');                        % MMP, 09/26/2026
  end
end

Pop = cdopvar(Cblk);
% % % END change MMP, 09/30/2026 (split)                                    % MMP, 09/30/2026

end


% % % BEGIN change MMP, 09/30/2026 (split): the local functions from here   % MMP, 09/30/2026
% to its END are the sections of the former primary, bodies unchanged.      % MMP, 09/30/2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = parse_inputs(dims,spaces,dom)                                  % MMP, 09/30/2026
% Spaces, registry, dims and domain, as the fields of S: M, nv, vars,       % MMP, 09/30/2026
% vars_int, vars_dum (the names of theta and s'), mask (M x nv, space k     % MMP, 09/30/2026
% has variable d), own{k} (its variable indices, a row), dims, dom.         % MMP, 09/30/2026

% % % BEGIN change MMP, 09/30/2026 (parser): spaces, registry, dims,        % MMP, 09/30/2026
% reserved names and domain are read by 'parse_copvar_spaces', the parser   % MMP, 09/30/2026
% this file shares with 'lpivar_cdopvar' and the _sop constructors. Deleted: % MMP, 09/30/2026
% the inline parse of the same forms (spaces; registry and mask; dims; the  % MMP, 09/30/2026
% reserved suffix check; domain; initial coding 09/21/2026, no later        % MMP, 09/30/2026
% markers), the parser's logic with six messages worded differently. A      % MMP, 09/30/2026
% 'cellstr' is still one space, so the single-space call matches            % MMP, 09/30/2026
% 'sopquadvar'. The registry comes back SORTED, which is not cosmetic: an   % MMP, 09/30/2026
% sopvar/sdopvar indexes its parameter cell over S3 in sorted order, and    % MMP, 09/30/2026
% with a sorted registry every space, and hence every pair's shared set,    % MMP, 09/30/2026
% is already sorted, so no cell permutation is needed at the end.           % MMP, 09/30/2026
[meta,spaces,sp_in] = parse_copvar_spaces(dims,spaces,dom);                 % MMP, 09/30/2026
% Pop is self-adjoint: the input spaces and dims must be the output ones.   % MMP, 09/30/2026
if ~isequal(sp_in,spaces) || ~isequal(meta.dim_in,meta.dim_out)             % MMP, 09/30/2026
    error("A self-adjoint operator has equal input and output spaces.")     % MMP, 09/30/2026
end                                                                         % MMP, 09/30/2026
M = numel(spaces);                                                          % MMP, 09/30/2026
vars = meta.vars;       nv = numel(vars);       mask = meta.space_out;      % MMP, 09/30/2026
dims = meta.dim_out;    dom = meta.dom;                                     % MMP, 09/30/2026

% Row orientation is forced on every direction list here and below. 'find'
% returns a COLUMN for a 1x1 input, so with a single registry variable an
% empty list comes back 0x1 rather than 1x0, and a 0x1 index silently turns
% a concatenation into an error and a one-row multi-index into a no-row one.
own = cell(1,M);
for k = 1:M
    own{k} = reshape(find(mask(k,:)),1,[]);
end

% % % Reserved names, as in 'sopquadvar': the integration variable theta and
% the input dummy s' are named by suffix, so the suffixes are reserved
% outright rather than only where they would collide.
% The check itself is 'check_reserved_names', called by the parser above.   % MMP, 09/30/2026
vars_int = strcat(vars,'_int');
vars_dum = strcat(vars,'_dum');
% % % END change MMP, 09/30/2026 (parser)                                   % MMP, 09/30/2026

S = struct('M',M,'nv',nv,'vars',{vars},'vars_int',{vars_int}, ...
           'vars_dum',{vars_dum},'mask',mask,'own',{own},'dims',dims,'dom',dom); % MMP, 09/30/2026

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [options,vartype,sep,incl] = parse_options(options,nv,M)           % MMP, 09/30/2026
% Options psatz (validated, as double), type, sep (1 x nv logical) and      % MMP, 09/30/2026
% include (1 x M cell), for nv registry variables and M spaces.             % MMP, 09/30/2026

% % % Options
% if nargin<6 || isempty(options)                                           % MMP, 09/30/2026 (was)
if isempty(options)             % the primary passes [] when omitted        % MMP, 09/30/2026
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if ~isfield(options,'psatz') || isempty(options.psatz)
    options.psatz = 0;
end
% Only 0 and 1 are implemented, and the value is checked rather than taken  % MMP, 09/27/2026 (was)
% 0, 1 and the face codes are implemented; the value is checked, not taken  % MMP, 09/27/2026
% for its truthiness. 'poslpivar_2d' also accepts 2, whose weight is a BALL
% (rds^2 - sum (s_d-c_d)^2) and not the box built below, so treating 2 as 1
% would hand a caller carrying a settings file across a different cone with
% no warning; 'settings2possopvar' drops psatz=2 terms for the same reason.
% if ~isscalar(options.psatz) || ~ismember(options.psatz,[0,1])             % MMP, 09/27/2026 (was)
%     error("'psatz' should be 0 or 1. 'poslpivar_2d' additionally accepts 2, "...
%           +"whose weight is a ball rather than the box this routine builds.") % MMP, 09/27/2026 (was)
% 2d+1 / 2d+2: one face of the box in registry direction d (see header).    % MMP, 09/27/2026
if ~isscalar(options.psatz) || ~ismember(options.psatz,[0,1,3:2*nv+2])      % MMP, 09/27/2026
    error("'psatz' should be 0, 1, or a face code 2d+1 / 2d+2 with d "...
          +"in 1..nv, here nv = "+num2str(nv)+"; 2 is "...
          +"poslpivar_2d's ball, not the box.")                             % MMP, 09/27/2026
end
% Double from here on: an integer class passes the check above, and then    % MMP, 09/27/2026
% floor((psatz-1)/2) rounds in integer arithmetic and picks the wrong face. % MMP, 09/27/2026
options.psatz = double(options.psatz);                                      % MMP, 09/27/2026
vartype = 'pos';
if isfield(options,'type') && ~isempty(options.type)
    vartype = char(options.type);
    if ~ismember(vartype,{'pos','sym'})
        error("'type' should be 'pos' or 'sym'.")
    end
end
if strcmp(vartype,'sym') && options.psatz
    error("'psatz' has no meaning for a 'sym' variable; use type 'pos'.")
end
if ~isfield(options,'sep') || isempty(options.sep)
    sep = false(1,nv);
else
    sep = logical(reshape(options.sep,1,[]));
    if isscalar(sep)
        sep = repmat(sep,1,nv);
    elseif numel(sep)~=nv
        error("'sep' should be a scalar or have one entry per registry variable.")
    end
end
if ~isfield(options,'include')
    incl = cell(1,M);
else
    incl = options.include;
    if ~iscell(incl) || iscellstr(incl)
        incl = {incl};              % one space given directly
    end
    incl = reshape(incl,1,[]);
    if numel(incl)~=M
        error("A cell 'include' should have one entry per space.")
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha_sp,nb,base,basis_list,bsp,bix] = basis_operators(S,incl)    % MMP, 09/30/2026
% The basis operators Z_alpha of Sec. 9.1, per space: alpha_sp{k} holds     % MMP, 09/30/2026
% one multi-index per row over space k's own variables, nb(k) rows; base,   % MMP, 09/30/2026
% bsp, bix map global basis c to (space, row) and back; basis_list is the   % MMP, 09/30/2026
% output of the same name.                                                  % MMP, 09/30/2026
M = S.M;    nv = S.nv;  own = S.own;    sep = S.sep;                        % MMP, 09/30/2026

% % % Basis operators of each space. The multi-index runs over the space's
% OWN variables only: a variable the space does not have carries no
% indicator, so there is nothing to choose in that direction.
alpha_sp = cell(1,M);   nb = zeros(1,M);
for k = 1:M
    nk = numel(own{k});
    vals = cell(1,nk);
    for t = 1:nk
        if sep(own{k}(t)),  vals{t} = [1,4];  else,  vals{t} = [1,2,3];  end
    end
%   A = alpha_grid(vals);                                                   % MMP, 09/30/2026 (was)
    % Shared enumeration, direction 1 fastest (Sec. 9.1 multi-indices): a   % MMP, 09/30/2026
    % linear 'include' index means the same here and in 'sopquadvar'.       % MMP, 09/30/2026
    A = multiindex_grid(vals);                                              % MMP, 09/30/2026
    A = alpha_select(A,incl{k},nk,k);
    if size(A,1)==0
        error("At least one basis operator must be included for space "...
              +num2str(k)+".")
    end
    if size(unique(A,'rows'),1)~=size(A,1)
        error("Multi-indices for space "+num2str(k)+" should be distinct.")
    end
    alpha_sp{k} = A;    nb(k) = size(A,1);
end
N = sum(nb);
base = cumsum([0,nb]);              % basis c of space k is base(k)+c globally

% Flat description of the basis list, and the reverse lookup used below.
basis_list = zeros(N,1+nv);
bsp = zeros(N,1);   bix = zeros(N,1);
for k = 1:M
    for i = 1:nb(k)
        c = base(k)+i;
        basis_list(c,1) = k;
        basis_list(c,1+own{k}) = alpha_sp{k}(i,:);
        bsp(c) = k;     bix(c) = i;
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deg_list = degree_list(deg,S,nb,bsp,bix)                           % MMP, 09/30/2026
% The degree specification of each basis operator: deg_list{c} is the       % MMP, 09/30/2026
% struct 'process_degrees_one' returns for global basis c, with fields      % MMP, 09/30/2026
%   int     1 x nv, cap on the degree in each integration variable theta_d; % MMP, 09/30/2026
%   mult   1 x nv, cap on the degree in each s_d (default: int);            % MMP, 09/30/2026
%   joint   cap on the total degree (default: no further cap);              % MMP, 09/30/2026
%   subset  caps per subset of [theta_1..theta_nv, s^k], or [];             % MMP, 09/30/2026
% the vocabulary of 'sopquadvar', over the registry. These cap the factor   % MMP, 09/30/2026
% Z^alpha(theta,s), not the kernel of Pop. Not the fields of the same       % MMP, 09/30/2026
% names in 'lpivar_cdopvar', which are kernel degrees per variable role:    % MMP, 09/30/2026
% there 'mult' is the left degree in a multiplier direction and 'int' the   % MMP, 09/30/2026
% [left right] degrees in an integral direction.                            % MMP, 09/30/2026
M = S.M;    nv = S.nv;  own = S.own;    N = sum(nb);                        % MMP, 09/30/2026

% % % Degrees, one struct per basis operator. 'int' and 'mult' are given over
% the registry; space k uses the 'mult' entries of its own variables.
if iscell(deg)
    if numel(deg)~=M
        error("A cell 'deg' should have one entry per space.")
    end
    deg_sp = reshape(deg,1,[]);
else
    deg_sp = repmat({deg},1,M);
end
deg_list = cell(1,N);
for c = 1:N
    k = bsp(c);     spec = deg_sp{k};
    if iscell(spec)
        if numel(spec)~=nb(k)
            error("A cell degree specification for space "+num2str(k)...
                  +" should have one entry per basis operator of that space.")
        end
        spec = spec{bix(c)};
    end
    d = process_degrees_one(spec,nv);
    if ~isempty(d.subset) && numel(own{k})~=nv
        % 'subset' is indexed over the variables the grid is built in, which
        % here is [theta_1..theta_nv, s^k]; that length is space-dependent,
        % so a registry-sized array would silently mean something else.
        error("A 'subset' degree array is only accepted for a space holding "...
              +"every registry variable; space "+num2str(k)+" does not.")
    end
    deg_list{c} = d;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Build the monomial bases Z^alpha(theta,s) and Z^alpha(theta,s')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Z1c,Z2c,mdim] = monomial_bases(S,alpha_sp,deg_list,bsp,bix)       % MMP, 09/30/2026
% Z1c{c}, Z2c{c}: the monomials of basis c in (theta,s) and (theta,s'),     % MMP, 09/30/2026
% the two arguments 'sosquadvar' takes; mdim(c) its space dimension.        % MMP, 09/30/2026
own = S.own;    vars = S.vars;  vars_int = S.vars_int;                      % MMP, 09/30/2026
vars_dum = S.vars_dum;  dims = S.dims;  N = numel(bsp);                     % MMP, 09/30/2026

% Exponents are stored over [theta_1..theta_nv, s^k], so Z1c{c} and Z2c{c}
% differ only in the name of the second group.
Z1c = cell(1,N);    Z2c = cell(1,N);    mdim = zeros(N,1);
for c = 1:N
    k = bsp(c);     ok = own{k};
    caps_int  = deg_list{c}.int;
    caps_mult = deg_list{c}.mult(ok);
    % A multiplier in direction d identifies s_d with theta_d, so any s_d
    % dependence of the basis would be redundant.
    caps_mult(alpha_sp{k}(bix(c),:)==1) = 0;
    Ei = build_exponent_grid([caps_int,caps_mult],deg_list{c}.joint, ...
                             deg_list{c}.subset);
    Ti = size(Ei,1);
    Z1c{c} = polynomial(speye(Ti),Ei,[vars_int(:);reshape(vars(ok),[],1)],[Ti,1]);
    Z2c{c} = polynomial(speye(Ti),Ei,[vars_int(:);reshape(vars_dum(ok),[],1)],[Ti,1]);
    mdim(c) = dims(k);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Declare Q and form the products N{c,e} = Z_c' Q_ce Z_e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The 'sosquadvar' call that declares Q is in the primary.                  % MMP, 09/30/2026

function [Zd,dmap,ndec] = common_dvars(Qcell)                               % MMP, 09/30/2026
% The sorted list Zd of every decision variable of Q, ndec = numel(Zd),     % MMP, 09/30/2026
% and dmap: name -> index into Zd.                                          % MMP, 09/30/2026

% Common decision variable basis, so that every block is expressed over the
% same Zd - a class invariant of 'cdopvar'.
% Pieces collected, concatenated once: growing Zd re-copied the whole list  % MMP, 09/26/2026
% for each of the N^2 blocks, O(N^2*q).                                     % MMP, 09/26/2026
% Zd = {};                                                                  % MMP, 09/26/2026 (was)
Zd = cell(numel(Qcell),1);                                                  % MMP, 09/26/2026
for q = 1:numel(Qcell)
%   Zd = [Zd; reshape(cellstr(string(Qcell{q})),[],1)];                     %#ok<AGROW> % MMP, 09/26/2026 (was)
    Zd{q} = reshape(cellstr(string(Qcell{q})),[],1);                        % MMP, 09/26/2026
end
Zd = vertcat(Zd{:});                                                        % MMP, 09/26/2026
Zd = unique(Zd);    ndec = numel(Zd);
dmap = containers.Map(Zd,num2cell(1:ndec));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gfun = psatz_weight(options,S)                                     % MMP, 09/30/2026
% The weight g(theta) of options.psatz (0: g = 1), see the help.            % MMP, 09/30/2026
nv = S.nv;  vars_int = S.vars_int;  dom = S.dom;                            % MMP, 09/30/2026

% Multiplier restricting positivity to the domain, evaluated at the
% integration variable as in 'poslpivar'.
gfun = polynomial(1);
% if options.psatz                                                          % MMP, 09/27/2026 (was)
if options.psatz==1                                                         % MMP, 09/27/2026
    for d = 1:nv
        thd = polynomial(vars_int(d));
        gfun = gfun*(thd-dom(d,1))*(dom(d,2)-thd);
    end
elseif options.psatz>=3                                                     % MMP, 09/27/2026
    % One face: odd code the lower face theta_d = ad, even the upper.       % MMP, 09/27/2026
    % Normalised by L_d, as in 'poslpivar_2d', so faces share a scale.      % MMP, 09/27/2026
    d = floor((options.psatz-1)/2);                                         % MMP, 09/27/2026
    thd = polynomial(vars_int(d));                                          % MMP, 09/27/2026
    if mod(options.psatz,2)==1                                              % MMP, 09/27/2026
        gfun = (thd-dom(d,1))/(dom(d,2)-dom(d,1));                          % MMP, 09/27/2026
    else                                                                    % MMP, 09/27/2026
        gfun = (dom(d,2)-thd)/(dom(d,2)-dom(d,1));                          % MMP, 09/27/2026
    end                                                                     % MMP, 09/27/2026
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Eliminate the integration variable, block pair by block pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [terms,nt] = block_terms(k,l,S,B)                                  % MMP, 09/30/2026
% The terms Z_i* Q_ij Z_j (Sec. 9.2) of block (k,l), i over the basis       % MMP, 09/30/2026
% operators of space k and j over those of space l, integration variable    % MMP, 09/30/2026
% eliminated: terms{1:nt}, one sdopvar per substitution pair of each        % MMP, 09/30/2026
% (i,j), all on the one list B.Zd. HOT: runs the whole elimination. It      % MMP, 09/30/2026
% must not assign into Zd, Nprod or prog: writing one element of a          % MMP, 09/30/2026
% shared q-length cell copies all of it (27.5 ms at q = 2.7e6, measured).   % MMP, 09/30/2026
nv = S.nv;  vars = S.vars;  vars_int = S.vars_int;  vars_dum = S.vars_dum;  % MMP, 09/30/2026
own = S.own;  dims = S.dims;  mask = S.mask;  dom = S.dom;  sep = S.sep;    % MMP, 09/30/2026
alpha_sp = B.alpha_sp;  nb = B.nb;  base = B.base;  Nprod = B.Nprod;        % MMP, 09/30/2026
gfun = B.gfun;  Zd = B.Zd;  dmap = B.dmap;  ndec = B.ndec;                  % MMP, 09/30/2026

% Adjoining a basis operator flips its lower and upper integrals; the
% full-domain integral is self-adjoint.
negmap = [1,3,2,4];

    ok = own{k};    ol = own{l};    mk = dims(k);   ml = dims(l);
%   nk = numel(ok); nlv = numel(ol);                                        % MMP, 09/30/2026 (was)
    nk = numel(ok);     % nlv went with the loops 'by_space' replaced       % MMP, 09/30/2026

    % Direction classes for this pair: 3 shared (S3), 1 output only (S2),
    % 2 input only (S1), 0 in neither space.
    cls = zeros(1,nv);
    cls( mask(k,:) &  mask(l,:)) = 3;
    cls( mask(k,:) & ~mask(l,:)) = 1;
    cls(~mask(k,:) &  mask(l,:)) = 2;
    D3 = reshape(find(cls==3),1,[]);        n3 = numel(D3);
    p3 = zeros(1,nv);       p3(D3) = 1:n3;      % position of d within D3

    vars_conv = struct('out',{[vars(ok),vars_int]},'in',{vars_dum(ol)});
    vars_io   = struct('out',{vars(ok)},'in',{vars(ol)});
    dom_io    = struct('out',dom(ok,:),'in',dom(ol,:));

    % A full-domain index expands into 2 substitutions per side, and one
    % arises only in a SEPARABLE direction, so the reachable bound is
    % 4^(separable shared directions), not 4^n3 -- which over-allocates an
    % empty cell array by 4^n3 per block pair with the default sep. This is
    % the bound 'sopquadvar' uses.
    terms = cell(nb(k)*nb(l)*4^sum(sep(D3)),1);     nt = 0;
    for i = 1:nb(k)
      for j = 1:nb(l)
        % % % Coefficients of the product, over the split (s,theta | s')
        Rij = gfun*Nprod{base(k)+i,base(l)+j};
        [Pblk,ZLb,ZRb] = dpvar2sdvar(Rij,vars_conv);

        Zs  = ZLb(1:nk);
        Zth = ZLb(nk+1:nk+nv);
        Zsp = ZRb;
        Ns  = prod([cellfun(@numel,Zs),1]);
        Nsp = prod([cellfun(@numel,Zsp),1]);
        NGa = prod([cellfun(@numel,Zth),1]);
        g1 = mk*Ns;     g2 = ml*Nsp;
        if Pblk.m~=g1*NGa || Pblk.n~=g2
            error("Internal error: unexpected coefficient dimensions.")
        end

        % Pack A and each row of B as sheets so one linear map serves them
        % all. The decision-variable axis is a pure BATCH axis, so the
        % elimination runs on this block's OWN rows and they are scattered
        % onto the global basis Zd once at the end; see the note at the
        % unpack loop below, and 'sopquadvar', which took the same change.
%       locB = dvar_rows(cellstr(string(Pblk.dvarname)),dmap);              % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
        % A cellstr needs no string round trip (0.6 s in 2-D Hinf).         % MMP, 09/26/2026
        dn = Pblk.dvarname;                                                 % MMP, 09/26/2026
        if ~iscellstr(dn),  dn = cellstr(string(dn));   end                 % MMP, 09/26/2026
        locB = dvar_rows(dn(:),dmap);                                       % MMP, 09/26/2026
        nloc = numel(locB);                                                 % MMP, 09/22/2026
        if size(Pblk.B,1)~=nloc                                             % MMP, 09/22/2026
            error("Internal error: B has %d rows for %d decision variables.",...% MMP, 09/22/2026
                  size(Pblk.B,1),nloc)                                      % MMP, 09/22/2026
        end                                                                 % MMP, 09/22/2026
%       Bblk = remap_dvars(Pblk.B,cellstr(string(Pblk.dvarname)),dmap,ndec);% MMP, 09/22/2026 (was)
%       Cpk  = pack_sheets([Pblk.A.';Bblk],Pblk.m,Pblk.n);                  % MMP, 09/22/2026 (was)
        Cpk  = pack_sheets([Pblk.A.';Pblk.B],Pblk.m,Pblk.n);                % MMP, 09/22/2026

        % % % Multi-indices of the pair, over the registry. The left factor
        % is the adjoint of basis i, so its indicators are flipped.
        aL = zeros(1,nv);   aL(ok) = negmap(alpha_sp{k}(i,:));
        aR = zeros(1,nv);   aR(ol) = alpha_sp{l}(j,:);
        oneidx = zeros(1,nv);
        oneidx(cls==1) = aL(cls==1);
        oneidx(cls==2) = aR(cls==2);

        % % % Directions that are not shared: one-sided or definite integrals
%       [Cpre,Zp,NL1,NR1] = int_onesided(Cpk,g1,g2,ndec,Zth,cls,oneidx,dom);% MMP, 09/22/2026 (was)
        [Cpre,Zp,NL1,NR1] = int_onesided(Cpk,g1,g2,nloc,Zth,cls,oneidx,dom);% MMP, 09/22/2026
        g1b = g1*NL1;       g2b = g2*NR1;

        % New exponents these produce, laid out per variable of each space so
        % they can be merged with the pre-existing bases.
%       ZpL = repmat({0},1,nk);                                             % MMP, 09/30/2026 (was)
%       for t = 1:nk                                                        % MMP, 09/30/2026 (was)
%           if cls(ok(t))==1,   ZpL{t} = Zp{ok(t)};     end                 % MMP, 09/30/2026 (was)
%       end                                                                 % MMP, 09/30/2026 (was)
%       ZpR = repmat({0},1,nlv);                                            % MMP, 09/30/2026 (was)
%       for t = 1:nlv                                                       % MMP, 09/30/2026 (was)
%           if cls(ol(t))==2,   ZpR{t} = Zp{ol(t)};     end                 % MMP, 09/30/2026 (was)
%       end                                                                 % MMP, 09/30/2026 (was)
        ZpL = by_space(Zp,ok,cls(ok)==1);                                   % MMP, 09/30/2026
        ZpR = by_space(Zp,ol,cls(ol)==2);                                   % MMP, 09/30/2026

        % % % Shared directions: the semiseparable integral. A full-domain
        % index expands into its 2/3 substitutions against the same Q.
        G = struct();
        G.C = Cpre;
        G.Z = Zth(D3);
        aL_all = expand_full(aL(D3));
        aR_all = expand_full(aR(D3));
        for iL = 1:size(aL_all,1)
        for iR = 1:size(aR_all,1)
%       [Cgam,ZL3,ZR3] = int_semisep(G,aL_all(iL,:),aR_all(iR,:),dom(D3,:));% MMP, 09/26/2026 (was)
        % Triplets, not matrices: each Cgam{q} would be g1b*NL3 x           % MMP, 09/26/2026
        % (1+nloc)*g2b*NR3, O(nloc) column pointers to write and then to    % MMP, 09/26/2026
        % scan in 'unpack_sheets' for a few thousand nonzeros.              % MMP, 09/26/2026
        [Cgam,ZL3,ZR3] = int_semisep(G,aL_all(iL,:),aR_all(iR,:), ...
                                     dom(D3,:),'triplets');                 % MMP, 09/26/2026

%       ZL3f = repmat({0},1,nk);                                            % MMP, 09/30/2026 (was)
%       for t = 1:nk                                                        % MMP, 09/30/2026 (was)
%           if cls(ok(t))==3,   ZL3f{t} = ZL3{p3(ok(t))};   end             % MMP, 09/30/2026 (was)
%       end                                                                 % MMP, 09/30/2026 (was)
%       ZR3f = repmat({0},1,nlv);                                           % MMP, 09/30/2026 (was)
%       for t = 1:nlv                                                       % MMP, 09/30/2026 (was)
%           if cls(ol(t))==3,   ZR3f{t} = ZR3{p3(ol(t))};   end             % MMP, 09/30/2026 (was)
%       end                                                                 % MMP, 09/30/2026 (was)
        ZL3f = by_space(ZL3,p3(ok),cls(ok)==3);                             % MMP, 09/30/2026
        ZR3f = by_space(ZR3,p3(ol),cls(ol)==3);                             % MMP, 09/30/2026

        % % % Merge the three groups of monomials on each side. The row index
        % is (m_k, Zs, ZpL, ZL3) and the column index (m_l, Zsp, ZpR, ZR3),
        % each group appended to the right of the previous one, so the maps
        % compose as kron(M1,I)*M2.
        [Z1m,ML1] = merge_monomial_product(Zs,ZpL);
        [ZLf,ML2] = merge_monomial_product(Z1m,ZL3f);
        NL3 = prod([cellfun(@numel,ZL3f),1]);
        ML  = kron(ML1,speye(NL3))*ML2;

        [W1m,MR1] = merge_monomial_product(Zsp,ZpR);
        [ZRf,MR2] = merge_monomial_product(W1m,ZR3f);
        NR3 = prod([cellfun(@numel,ZR3f),1]);
        MR  = kron(MR1,speye(NR3))*MR2;

        Lmat = kron(speye(mk),ML).';
        Rmat = kron(speye(ml),MR);

        % % % Parameter cell of the term over the gamma cells of Sec. 6.1   % MMP, 09/30/2026
        params = gamma_params(Cgam,g1b,NL3,g2b,NR3,nloc,Lmat,Rmat, ...
                              locB,ndec,n3);                                % MMP, 09/30/2026

        nt = nt+1;
        terms{nt} = sdopvar(params,vars_io,Zd,ZLf,ZRf,dom_io,[mk,ml]);
        end
        end
      end
    end
    if nt==0
        error("Internal error: block (%d,%d) collected no terms.",k,l)
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = gamma_params(Cgam,g1b,NL3,g2b,NR3,nloc,Lmat,Rmat, ...
                               locB,ndec,n3)                                % MMP, 09/30/2026
% The sdopvar parameter cell of one term: params.A{q}, params.B{q} (the     % MMP, 09/30/2026
% constant and decision parts of gamma cell q, B over the global list of    % MMP, 09/30/2026
% ndec variables), as a 3 x ... x 3 cell over the n3 shared variables.      % MMP, 09/30/2026
% Cgam{q}: triplets from 'int_semisep' on the block's own nloc sheets;      % MMP, 09/30/2026
% locB: their rows in the global list; Lmat, Rmat: the basis merge maps.    % MMP, 09/30/2026

        % Nothing between 'pack_sheets' and here mixes sheets -- each stays a
        % contiguous block of columns, as 'unpack_sheets' own header records --
        % so the block's own nloc sheets carry the elimination and the rows are
        % scattered onto the global basis once, here, from triplets. Widening
        % to ndec sheets beforehand padded every intermediate: measured in
        % 'sopquadvar' at 134 to 2119 times the width the block needs.
        params = struct();
        params.A = cell(numel(Cgam),1);
        params.B = cell(numel(Cgam),1);
        % Lmat and Rmat are those of this (iL,iR) pair for every gamma cell, % MMP, 09/26/2026
        % so their Kronecker product is formed once here, not once per cell. % MMP, 09/26/2026
%       KLR = kron(Rmat.',Lmat).';                                          % MMP, 09/26/2026 (was)
        % Superseded the same day: 'unpack_sheets' now applies Lmat, Rmat and
        % the scatter onto Zd as one index remap of the triplets of Cgam{q}, so
        % no Kronecker product is formed at all (see its header).
        for q = 1:numel(Cgam)
%           [Aq,Bq] = unpack_sheets(Cgam{q},g1b*NL3,g2b*NR3,ndec);          % MMP, 09/22/2026 (was)
%           [Aq,Bq] = unpack_sheets(Cgam{q},g1b*NL3,g2b*NR3,nloc);          % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
%           [params.A{q},Bout] = lr_multiply(Lmat,Aq,Bq,Rmat);              % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
%           [params.A{q},Bout] = lr_multiply(Lmat,Aq,Bq,Rmat,KLR);          % MMP, 09/26/2026 (was)
%           [params.A{q},params.B{q}] = lr_multiply(Lmat,Aq,Bq,Rmat);       % MMP, 09/22/2026 (was)
            % 'find' returns ROWS for a single-row input, which Bout is when
            % the block carries one decision variable; forced to columns, as
            % elsewhere in this file.
%           [bi,bj,bv] = find(Bout);                                        % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
%           params.B{q} = sparse(colv(locB(bi)),colv(bj),colv(bv), ...
%                                ndec,size(Bout,2));                        % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
            % All-zero gamma cell: write its zeros, skip the unpack. The    % MMP, 09/29/2026
            % Sec. 6.1 key table admits at most 2 gammas per direction for  % MMP, 09/29/2026
            % one (beta,alpha), so at most 2^n3 of the 3^n3 cells are       % MMP, 09/29/2026
            % nonzero (91% of w3 cells empty). 'triplets' stores no zero    % MMP, 09/29/2026
            % v, so an empty v is an all-zero cell, and these are the       % MMP, 09/29/2026
            % zeros 'unpack_sheets' returns for it: same size and class.    % MMP, 09/29/2026
            % The skip 'sopquadvar' has had since 09/22/2026.               % MMP, 09/29/2026
            if isempty(Cgam{q}.v)                                           % MMP, 09/29/2026
                nAB = size(Lmat,1)*size(Rmat,2);                            % MMP, 09/29/2026
                params.A{q} = sparse(nAB,1);                                % MMP, 09/29/2026
                params.B{q} = sparse(ndec,nAB);                             % MMP, 09/29/2026
                continue                                                    % MMP, 09/29/2026
            end                                                             % MMP, 09/29/2026
            [params.A{q},params.B{q}] = unpack_sheets(Cgam{q},g1b*NL3, ...
                g2b*NR3,nloc,Lmat,Rmat,locB,ndec);                          % MMP, 09/26/2026
        end
        % The parameter cell is indexed over the pair's shared variables. The
        % registry is sorted, so D3 is already in the sorted order the class
        % reads a cell subscript in and no permutation is needed.
        params.A = reshape(params.A,[3*ones(1,n3),1,1]);
        params.B = reshape(params.B,[3*ones(1,n3),1,1]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Z = by_space(src,idx,sel)                                          % MMP, 09/30/2026
% Exponent lists laid out per variable of one space: Z{t} = src{idx(t)}     % MMP, 09/30/2026
% where sel(t), and 0 (degree 0 only) elsewhere. Replaces four identical    % MMP, 09/30/2026
% scatter loops of 'block_terms' (initial coding 09/21/2026).               % MMP, 09/30/2026

Z = repmat({0},1,numel(idx));                                               % MMP, 09/30/2026
Z(sel) = src(idx(sel));                                                     % MMP, 09/30/2026

end
% % % END change MMP, 09/30/2026 (split)                                    % MMP, 09/30/2026


% MMP, 09/30/2026: deleted the local function 'alpha_grid' (initial coding  % MMP, 09/30/2026
% 09/21/2026, no later markers), now 'multiindex_grid' (misc/conventions),  % MMP, 09/30/2026
% whose body is alpha_grid's.                                               % MMP, 09/30/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = alpha_select(A,incl,nd,k)
% Apply one space's 'include' specification, in any of the three formats
% 'sopquadvar' accepts. Validating multi-index rows against the generated
% list rather than a numeric range is what makes a separable direction, which
% admits only 1 and 4, reject 2 and 3.

if isempty(incl)
    return
end
if islogical(incl)
    if numel(incl)~=size(A,1)
        error("A logical 'include' for space "+num2str(k)+" has one entry "...
              +"per basis operator of that space.")
    end
    A = A(incl(:),:);
elseif nd>0 && size(incl,2)==nd
    if ~all(ismember(incl,A,'rows'))
        error("'include' for space "+num2str(k)+" contains an inadmissible "...
              +"multi-index.")
    end
    A = incl;
else
    if any(incl(:)<1) || any(incl(:)>size(A,1)) || any(incl(:)~=round(incl(:)))
        error("Linear indices in 'include' for space "+num2str(k)...
              +" are out of range.")
    end
    A = A(incl(:),:);
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C2,Zp,NL1,NR1] = int_onesided(C,g1,g2,nsheet,Zth,cls,oneidx,dom)
% Integrate out the theta directions that are NOT shared by the block pair.
%
% Sharing is what makes the semiseparable integral necessary: with a variable
% on both sides the result carries an indicator in s-s' and splits over
% gamma. With the variable on one side only, or on neither, the integral is
% an antiderivative of a monomial evaluated at one limit,
%
%   cls(d)==1   int I_beta(s_d-theta_d) theta_d^p dtheta_d,  beta = oneidx(d)
%   cls(d)==2   int I_alpha(theta_d-s'_d) theta_d^p dtheta_d
%   cls(d)==0   int_a^b theta_d^p dtheta_d
%
% and the outcome is a polynomial in s_d, in s'_d, or a constant. Doing these
% first leaves 'int_semisep' only the shared directions.
%
% INPUTS
% - C:      (g1*NG) x ((1+nsheet)*g2) packed coefficient matrix, rows indexed
%           (g1, theta) with theta fastest and direction 1 slowest within it,
%           columns (sheet, g2);
% - Zth:    1 x nv cell of theta exponent column vectors;
% - cls:    1 x nv direction classes, 3 shared and left alone;
% - oneidx: 1 x nv multi-index entries, read only where cls is 1 or 2;
%
% OUTPUTS
% - C2:     (g1*NL1*NG3) x ((1+nsheet)*g2*NR1), rows now indexed
%           (g1, new s monomials, remaining theta) and columns
%           (sheet, g2, new s' monomials), which is the layout
%           'int_semisep' consumes with g1*NL1 and g2*NR1 as its dimensions;
% - Zp:     1 x nv cell, Zp{d} the new exponents produced in direction d and
%           [] where none are;
% - NL1,NR1: sizes of the two new monomial groups.

nv = numel(cls);
W = cell(1,nv);     nout = ones(1,nv);      Zp = cell(1,nv);
for d = 1:nv
    P = reshape(Zth{d},[],1);
    switch cls(d)
        case 3
            W{d} = speye(numel(P));
        case 0
            [~,W{d}] = definite_int(P,dom(d,:));
        otherwise
            [Zp{d},W{d}] = onesided_int(P,cls(d),oneidx(d),dom(d,:));
    end
    nout(d) = size(W{d},1);
end

NG  = prod([cellfun(@numel,Zth),1]);
dL  = reshape(find(cls==1),1,[]);
dR  = reshape(find(cls==2),1,[]);
d3  = reshape(find(cls==3),1,[]);
NL1 = prod([nout(dL),1]);
NR1 = prod([nout(dR),1]);
NG3 = prod([nout(d3),1]);

% One transfer table for all directions at once. Direction 1 is the slowest
% index on both sides, which is what makes the Kronecker product in this
% order agree with the row layout of C.
T = 1;
for d = 1:nv
    T = kron(T,W{d});
end
% 'find' returns ROW vectors for a single-row input, which T is whenever
% every direction collapses, so the orientation is forced rather than
% assumed: a row 'tc' would be read by 'accumarray' below as one
% n-dimensional subscript.
[tr,tc,tw] = find(T);
tr = tr(:);     tc = tc(:);     tw = tw(:);

% Split the transfer's output index into the three groups, each with the
% lowest-numbered direction slowest so that it matches a per-variable
% monomial cell built in registry order.
% The split is the kron order of 'kron_split' (misc/conventions), kept      % MMP, 09/30/2026
% inline: the helper calls made this function 5 us (2%) slower per call,   % MMP, 09/30/2026
% measured on 51 captured w3 calls, and it runs 3248 times per w3 build.    % MMP, 09/30/2026
ns = numel(tr);
sub = zeros(ns,nv);
rem = tr-1;
for d = nv:-1:1
    sub(:,d) = mod(rem,nout(d));
    rem = floor(rem/nout(d));
end
iL = zeros(ns,1);   for t = 1:numel(dL),  iL = iL*nout(dL(t))+sub(:,dL(t));  end
iR = zeros(ns,1);   for t = 1:numel(dR),  iR = iR*nout(dR(t))+sub(:,dR(t));  end
i3 = zeros(ns,1);   for t = 1:numel(d3),  i3 = i3*nout(d3(t))+sub(:,d3(t));  end

% Apply it by remapping indices. A row of C moves partly to a column, since a
% direction that only the input space has produces an s' monomial, so this is
% not a left or right multiplication and is done on the triplets. The cost is
% the number of nonzeros times the fan-out of the table; nothing densifies
% along the decision variable axis, which is the sheet index of the columns.
[ir,ic,v] = find(C);
ir = ir(:);     ic = ic(:);     v = v(:);       % see the note on T above
arow = floor((ir-1)/NG);            % 0-based (matrix row, s monomial)
ith  = mod(ir-1,NG)+1;
sh   = floor((ic-1)/g2);            % 0-based sheet
jc   = mod(ic-1,g2);                % 0-based (matrix column, s' monomial)

cnt = accumarray(tc,1,[NG,1]);
ptr = [0;cumsum(cnt)];
[~,ord] = sort(tc);
iLs = iL(ord);  iRs = iR(ord);  i3s = i3(ord);  tws = tw(ord);

cc  = cnt(ith);
tot = sum(cc);
if tot==0
    C2 = sparse(g1*NL1*NG3,(1+nsheet)*g2*NR1);
    return
end
% Every 'repelem' result is forced to a column. 'repelem' of a 1 x 1 first
% argument returns a ROW, and 'pos' below is a column, so when C has exactly
% ONE nonzero the subtraction implicitly expands into a tot x tot matrix
% instead of a tot x 1 vector. Every entry then gets emitted 'cc' times at
% the SAME (row,column), 'sparse' sums the duplicates, and the block comes
% out scaled by cc with no error anywhere.
%
% MEASURED: for spaces {s1,s2} x {s2} at degree 0 with one basis each,
% T = [2;-1] and a single nonzero in C gave sel = [1 2 1 2] and a kernel of
% 2*(b1-s1) against the correct (b1-s1). It bites only when nnz(C)==1 AND a
% non-shared direction has fan-out above 1, which is why it is invisible at
% degree 1 or more (C then has several nonzeros), invisible for a multiplier
% direction (fan-out 1), and was invisible on a domain starting at 0, where
% the alpha=2 weight -a^(p+1)/(p+1) is zero and drops out leaving fan-out 1.
% The span comparisons against poslpivar cannot see it either: scaling a
% generator does not change the span.
off = cumsum([0;cc(1:end-1)]);
pos = (1:tot).' - colv(repelem(off,cc));        % 1..cc(u) within nonzero u
sel = colv(repelem(ptr(ith),cc)) + pos;

newrow = (colv(repelem(arow,cc))*NL1 + iLs(sel))*NG3 + i3s(sel) + 1;
newcol = colv(repelem(sh,cc))*(g2*NR1) + colv(repelem(jc,cc))*NR1 + iRs(sel) + 1;
newval = colv(repelem(v,cc)).*tws(sel);

C2 = sparse(newrow,newcol,newval,g1*NL1*NG3,(1+nsheet)*g2*NR1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = colv(x)
% Force a column, so that an arithmetic combination of two such vectors
% cannot silently become a matrix by implicit expansion.

x = reshape(x,[],1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zo,W] = onesided_int(P,side,idx,lim)
% Integral of theta^p against a single indicator over the domain, returning
% the exponents Zo of the result in the surviving spatial variable and the
% map W from the theta exponents P to them.
%
%   side 1   the variable is the OUTPUT s_d, indicator I_idx(s_d-theta_d)
%   side 2   the variable is the INPUT s'_d, indicator I_idx(theta_d-s'_d)
%
% so the limit carrying the variable depends on both, since the indicator's
% argument is negated between the two sides:
%
%   side 1, idx 2   I_1(s-th)=1 for th<=s    ->  int_a^s
%   side 1, idx 3   I_-1(s-th)=1 for th>=s   ->  int_s^b
%   side 2, idx 2   I_1(th-s')=1 for th>=s'  ->  int_{s'}^b
%   side 2, idx 3   I_-1(th-s')=1 for th<=s' ->  int_a^{s'}

if idx==1
    % delta identifies theta_d with the spatial variable: substitution.
    Zo = reshape(P,[],1);
    W = speye(numel(P));
    return
end
if idx==4
    % A full-domain integral has constant limits, so nothing survives in the
    % spatial variable.
    [Zo,W] = definite_int(P,lim);
    return
end

a = lim(1);     b = lim(2);     P = reshape(P,[],1);    np = numel(P);
var_is_upper = (side==1 && idx==2) || (side==2 && idx==3);

Zo = unique([0;P+1]);           % P>=0, so Zo(1)==0 and the constant row is 1
[~,posv] = ismember(P+1,Zo);
w0 = zeros(np,1);   wv = zeros(np,1);
for t = 1:np
    p = P(t);
    if var_is_upper
        % int_a^x th^p dth = (x^(p+1) - a^(p+1))/(p+1)
        w0(t) = -a^(p+1)/(p+1);     wv(t) =  1/(p+1);
    else
        % int_x^b th^p dth = (b^(p+1) - x^(p+1))/(p+1)
        w0(t) =  b^(p+1)/(p+1);     wv(t) = -1/(p+1);
    end
end
W = sparse([ones(np,1);posv],[(1:np).';(1:np).'],[w0;wv],numel(Zo),np);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zo,W] = definite_int(P,lim)
% int_a^b theta^p dtheta for each exponent in P: a single constant row.

a = lim(1);     b = lim(2);     P = reshape(P,[],1);    np = numel(P);
Zo = 0;
w = zeros(1,np);
for t = 1:np
    p = P(t);
    w(t) = (b^(p+1)-a^(p+1))/(p+1);
end
W = sparse(ones(1,np),1:np,w,1,np);

end
