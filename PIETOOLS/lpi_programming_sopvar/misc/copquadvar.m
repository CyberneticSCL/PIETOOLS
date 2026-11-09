function [prog,Pop,Qcell,basis_list] = copquadvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,BASIS_LIST] = COPQUADVAR(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares a self-adjoint 'cdopvar' decision operator Pop on the product of
% M spaces
%
%       X = X_1 x ... x X_M,   X_k = L_2^{m_k}[s^k]  (R^{m_k} if s^k = {}),
%
% as the quadratic form (sopvar document, Sec. 8.4 and 9)
%
%       Pop = sum_{c,e} (Z_c)* (g Q_{ce}) (Z_e),
%
% c and e over the basis operators of all spaces, with one Gram matrix
% Q = [Q_{ce}] over all of them, symmetric, and positive semidefinite for
% type 'pos'; block (k,l) of Pop is the sum over the c of space k and the e
% of space l. Basis operator c of space k carries a multi-index alpha over
% the variables of s^k and maps X_k to L_2^{m_k}[theta], theta =
% (theta_1,...,theta_nv) the variables of the registry, the sorted union of
% the s^k:
%
%       (Z_c x)(theta) = int I_alpha(theta-s) (I_{m_k} kron Z^c(theta,s))
%                            x(s) ds,
%
%       I_alpha(r) = prod_{d in s^k} I_{alpha_d}(r_d),
%       I_0(r) = delta(r),  I_1(r) = 1 for r >= 0,  I_-1(r) = 1 for r <= 0,
%
% the integral over the domain of s^k (none on R^{m_k}), I_1 and I_-1 zero
% elsewhere. A multi-index stores 1, 2, 3 for the entries 0 (multiplier), 1
% (lower integral) and -1 (upper integral), and 4 for the full-domain
% integral I_1 + I_-1 = 1. Z^c(theta,s) is the vector of the monomials in
% theta and s^k within the degree caps of c (degree 0 in s_d where alpha_d
% is a multiplier), and g(theta) the weight of OPTIONS.psatz, a multiplier
% on L_2[theta]. Integrating theta out, a direction d of s^k and s^l gives
% a semiseparable kernel in (s_d,s'_d), one of s^k only a polynomial in
% s_d, one of s^l only a polynomial in s'_d, and one of neither a definite
% integral; each block of Pop is an 'sdopvar' object, affine in the entries
% of Q.
%
% INPUT
% - prog:    'struct', the LPI program to which the decision variables of Q
%            are added;
% - dims:    'double', M x 1 array of positive integers, dims(k) = m_k, or
%            a scalar for every space; or a 'struct' with field out (that
%            array) and an optional field in equal to out;
% - spaces:  1 x M 'cell', spaces{k} the variables s^k of space k: a
%            'cellstr' of distinct names, a 'char' name, a 'polynomial'
%            object (its variables), or empty for R^{m_k}. A 'cellstr',
%            'char' or 'polynomial' SPACES is one space, an empty SPACES
%            the space R^{m_1}. Or a 'struct' with field out (that cell)
%            and an optional field in equal to out. Names ending in '_int'
%            or '_dum' are an error;
% - dom:     'double', nv x 2 array, dom(d,:) = [a_d,b_d], a_d < b_d, the
%            domain of registry variable d; or a 1 x 2 row for every
%            variable; or a 'struct' with fields vars (a 'cellstr' of the
%            names) and dom (one row per name);
% - deg:     the degree caps of the monomials Z^c(theta,s), over the
%            variables [theta_1..theta_nv, s^k]: a 'double' scalar d, the
%            cap of every theta_d and every s_d; or a 'struct' with the
%            optional fields
%              int     scalar or 1 x nv, the cap in each theta_d (default
%                      1);
%              mult    scalar or 1 x nv over the registry, the cap in each
%                      s_d of s^k (default int);
%              joint   scalar, the cap of the total degree (default none);
%              subset  2^(2*nv) array, the cap of the total degree in each
%                      subset of [theta_1..theta_nv, s_1..s_nv], the subset
%                      with members b at index 1+sum(2.^(b-1)); accepted
%                      only for a space holding every registry variable
%                      (default none);
%            or a 1 x M 'cell', one scalar or 'struct' per space, or for a
%            space a 'cell' with one per basis operator of that space, in
%            the order of its rows of BASIS_LIST;
% - options: (optional) 'struct' with the optional fields
%              type     'char', 'pos' (default) for Q >= 0, so that
%                       Pop = Pop* >= 0, or 'sym' for Q symmetric;
%              psatz    'double', the weight g: 0 (default) g = 1; 1 g =
%                       prod_d (theta_d-a_d)(b_d-theta_d); 2d+1 or 2d+2,
%                       d = 1..nv over the registry, g = (theta_d-a_d)/
%                       (b_d-a_d) or (b_d-theta_d)/(b_d-a_d). Any other
%                       value is an error, and a nonzero one with type
%                       'sym';
%              sep      'logical' scalar or 1 x nv over the registry
%                       (default false): where sep(d), alpha_d takes the
%                       values 1 and 4, elsewhere 1, 2 and 3;
%              include  1 x M 'cell', include{k} the basis operators of
%                       space k among all its multi-indices (one column
%                       per variable of s^k in registry order, enumerated
%                       with the first column fastest): an array of rows
%                       of multi-indices, a 'logical' mask, or linear
%                       indices into that enumeration; [] (default) for
%                       all; for M = 1 the entry itself. Each space keeps
%                       at least one, all distinct.
%
% OUTPUT
% - prog:       'struct', PROG with the decision variables of Q and, for
%               type 'pos', the constraint Q >= 0;
% - Pop:        'cdopvar' object, M x M, on X: block (k,l) an 'sdopvar'
%               object from X_l to X_k, every block on the one list of the
%               decision variables of Q (sorted names); Pop = Pop*, and
%               Pop >= 0 for type 'pos';
% - Qcell:      N x N 'cell', N the number of basis operators of all
%               spaces: Qcell{c,e} the 'cellstr' array of the names of the
%               entries of Q_{ce}, rows (component, monomial of Z^c) and
%               columns (component, monomial of Z^e), monomial fastest;
% - basis_list: N x (1+nv) 'double', row c for basis operator c: its space
%               index k in column 1 and its multi-index entry alpha_d in
%               column 1+d (0 where d is not a variable of s^k); the rows
%               ordered by space, then as included.
%
% See also SOPQUADVAR, POSCOPVAR, POSSOPVAR, SOSQUADVAR, INT_SEMISEP,
% CDOPVAR, LPI_EQ_CDOPVAR.
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
% MMP, 09/21/2026: Every 'repelem' result in 'int_onesided' forced to a
%                  column: 'repelem' of a scalar returns a row, so a block
%                  with one nonzero was scaled by the one-sided fan-out
%                  (reachable at monomial degree 0 only; invisible to a span
%                  comparison). options.psatz validated by value (2 built
%                  the box where 'poslpivar_2d' builds a ball). The per-pair
%                  term preallocation bounded by the separable shared
%                  directions, not all shared ones (4^n3 over-allocation).
% MMP, 09/22/2026: The scatter onto the global decision list deferred to
%                  after the elimination, as 'sopquadvar' the same day:
%                  'int_onesided' receives nloc, not ndec, as its sheet
%                  count (5.8x at degree 2, ndec = 1.7e5; coefficients
%                  unchanged, 'test_poscopvar').
% MMP, 09/25/2026: Container classes renamed mopvar -> copvar, mdopvar ->
%                  cdopvar; here mopquadvar -> copquadvar, test_posmopvar ->
%                  test_poscopvar. File was 'mopquadvar.m'.
% MMP, 09/26/2026: The per-gamma unpack, 'lr_multiply', find and scatter
%                  onto Zd are one 'unpack_sheets' call applying Lmat and
%                  Rmat as an index remap of the triplets of Cgam{q} (no
%                  Kronecker product formed); int_semisep called with layout
%                  'triplets', so the (1+nloc)*g2b*NR3-column matrix is
%                  never built; Zd concatenated once; plus_batch(...,
%                  'shared_Zd') skips its O(q) per-operand list comparison.
%                  Programs bit-identical; 2-D H-infinity build 55 -> 40 s.
% MMP, 09/27/2026: options.psatz = 2d+1 / 2d+2 weights the Gram form by the
%                  normalised linear generator of one face of the box, d
%                  over the registry (the N-D extension of poslpivar_2d
%                  psatz 3-6): the product weight vanishes on every face,
%                  so it cannot certify near one (2-D heatNd: the product
%                  certifies nothing, the 4 faces reach 0.9875 kappa*).
% MMP, 09/29/2026: An all-zero gamma cell skipped before 'unpack_sheets' and
%                  its zeros written directly (at most 2^n3 of the 3^n3
%                  cells can be nonzero, Sec. 6.1); programs bit-identical.
%                  An empty B stored as [] is rejected by 'plus_batch'.
% MMP, 09/30/2026: Shared helpers (sopvar/misc/conventions): spaces,
%                  registry, dims, reserved names and domain read by
%                  'parse_copvar_spaces' (the self-adjoint check stays
%                  here); 'multiindex_grid' replaces the local 'alpha_grid'.
%                  Same values on every input the old parse accepted; six
%                  messages reworded, a struct array 'spaces' refused, ''
%                  one R^q space, struct dims and [] spaces accepted.
% MMP, 09/30/2026: The primary split into local functions, one per phase,
%                  bodies moved unchanged: parse_inputs, parse_options,
%                  basis_operators, degree_list, monomial_bases,
%                  common_dvars, psatz_weight, block_terms (with
%                  gamma_params and by_space). Two moved lines change,
%                  marked (was): the nargin test on options, and the unused
%                  nlv. Programs bit-identical, time unchanged. The help and
%                  degree_list state what each 'deg' field means here and
%                  that 'lpivar_cdopvar' gives 'int' and 'mult' other
%                  meanings.
% MMP, 10/01/2026: 'sopquadvar' calls this routine (one space) instead of
%                  its own pair loop; no code change here.
% MMP, 10/28/2026: Header to the header protocol, step comments in the body
%                  and headers of the local functions; code unchanged.
%                  Deleted: the comparisons with sopquadvar, poslpivar,
%                  poslpivar_2d and lpivar_cdopvar, the history of the single
%                  implementation and of the shared parser, a measurement
%                  and the rationale for the one Gram and the order of the
%                  integrations.
% MMP, 11/08/2026: Comment pass (maintainer, 11/08/2026): log condensed and
%                  in date order, a (was) comment line and a note on (was)
%                  code removed, the notes in gamma_params and int_onesided
%                  shortened, See also without the retired poslpivar; code
%                  unchanged.
% MMP, 11/09/2026: Inlining pass (maintainer, 11/09/2026): the local
%                  COMMON_DVARS (09/30/2026, one caller) inlined where B.Zd,
%                  B.dmap and B.ndec are set, its header as the comment of
%                  the block; its 09/26/2026 lines keep their markers.


% % % BEGIN change MMP, 09/30/2026 (split): the primary runs the            % MMP, 09/30/2026
% construction of Sec. 9 phase by phase, one call each: the inputs; the     % MMP, 09/30/2026
% basis operators Z_alpha (Sec. 9.1) and their degree caps; their           % MMP, 09/30/2026
% monomial bases; one Gram Q through 'sosquadvar' (Sec. 9.3); the pairs     % MMP, 09/30/2026
% Z_i* Q_ij Z_j (Sec. 9.2), block by block. Each phase was a section of     % MMP, 09/30/2026
% this function and is now the local function below named after it, its     % MMP, 09/30/2026
% lines moved unchanged with their markers. New, and marked: the lines      % MMP, 09/30/2026
% passing variables in and out, and in 'block_terms' the calls to           % MMP, 09/30/2026
% 'gamma_params' and 'by_space'.                                            % MMP, 09/30/2026

% Five inputs required; OPTIONS [] when omitted
if nargin<5
    error("Not enough input arguments.")
end
if nargin<6,    options = [];   end                                         % MMP, 09/30/2026

% Spaces and registry; options; basis operators, their degree caps and
% their monomial vectors
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
% [B.Zd,B.dmap,B.ndec] = common_dvars(Qcell);                               % MMP, 11/09/2026 (was); MMP, 09/30/2026
% The decision variables of the Gram blocks Qcell (the former local         % MMP, 11/09/2026
% COMMON_DVARS): B.Zd the sorted list of every name, without repeats, the   % MMP, 11/09/2026
% common basis every block is expressed over (a class invariant of          % MMP, 11/09/2026
% 'cdopvar'); B.ndec = numel(B.Zd); B.dmap: name -> index into B.Zd.        % MMP, 11/09/2026
% Pieces collected, concatenated once: growing the list re-copied it for    % MMP, 11/09/2026; MMP, 09/26/2026
% each of the N^2 blocks, O(N^2*q).                                         % MMP, 11/09/2026; MMP, 09/26/2026
% Zd = {};                                                                  % MMP, 11/09/2026; MMP, 09/26/2026 (was)
Zd = cell(numel(Qcell),1);                                                  % MMP, 11/09/2026; MMP, 09/26/2026
for q = 1:numel(Qcell)                                                      % MMP, 11/09/2026
%   Zd = [Zd; reshape(cellstr(string(Qcell{q})),[],1)];                     %#ok<AGROW> % MMP, 11/09/2026; MMP, 09/26/2026 (was)
    Zd{q} = reshape(cellstr(string(Qcell{q})),[],1);                        % MMP, 11/09/2026; MMP, 09/26/2026
end                                                                         % MMP, 11/09/2026
Zd = vertcat(Zd{:});                                                        % MMP, 11/09/2026; MMP, 09/26/2026
B.Zd = unique(Zd);    B.ndec = numel(B.Zd);                                 % MMP, 11/09/2026
B.dmap = containers.Map(B.Zd,num2cell(1:B.ndec));                           % MMP, 11/09/2026
B.gfun = psatz_weight(options,S);                                           % MMP, 09/30/2026

% Block (k,l): the sum of the terms of all pairs of basis operators of
% spaces k and l
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

% The container of the M x M blocks
Pop = cdopvar(Cblk);
% % % END change MMP, 09/30/2026 (split)                                    % MMP, 09/30/2026

end


% % % BEGIN change MMP, 09/30/2026 (split): the local functions from here   % MMP, 09/30/2026
% to its END are the sections of the former primary, bodies unchanged.      % MMP, 09/30/2026

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Process the inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = parse_inputs(dims,spaces,dom)                                  % MMP, 09/30/2026
% S = PARSE_INPUTS(DIMS,SPACES,DOM) reads the spaces, dimensions and domain
% of COPQUADVAR, after checking that the input spaces and dimensions are
% the output ones.
% Spaces, registry, dims and domain, as the fields of S: M, nv, vars,       % MMP, 09/30/2026
% vars_int, vars_dum (the names of theta and s'), mask (M x nv, space k     % MMP, 09/30/2026
% has variable d), own{k} (its variable indices, a row), dims, dom.         % MMP, 09/30/2026
% INPUT  dims, spaces, dom: the inputs of COPQUADVAR of these names.
% OUTPUT S: 'struct' with the fields above: M, nv 'double' (the numbers of
%           spaces and of registry variables); vars (sorted), vars_int,
%           vars_dum 1 x nv 'cellstr'; mask M x nv 'logical'; own 1 x M
%           'cell' of 'double' rows; dims M x 1 'double'; dom nv x 2
%           'double'.

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
% Number of spaces; the sorted registry, the variable masks, dims, domain
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
% [OPTIONS,VARTYPE,SEP,INCL] = PARSE_OPTIONS(OPTIONS,NV,M) checks the
% options of COPQUADVAR and returns them with their defaults:
% Options psatz (validated, as double), type, sep (1 x nv logical) and      % MMP, 09/30/2026
% include (1 x M cell), for nv registry variables and M spaces.             % MMP, 09/30/2026
% INPUT  options: 'struct', the OPTIONS of COPQUADVAR, or [] when omitted;
%        nv, M:   'double', the numbers of registry variables and spaces.
% OUTPUT options: OPTIONS with field psatz, a 'double' (0 when absent);
%        vartype: 'char', 'pos' or 'sym';
%        sep:     1 x nv 'logical';
%        incl:    1 x M 'cell', the include entry of each space ([]: all).

% % % Options
% if nargin<6 || isempty(options)                                           % MMP, 09/30/2026 (was)
if isempty(options)             % the primary passes [] when omitted        % MMP, 09/30/2026
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
% psatz: 0 when absent
if ~isfield(options,'psatz') || isempty(options.psatz)
    options.psatz = 0;
end
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
% type: 'pos' when absent, or 'sym'; a psatz weight only for 'pos'
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
% sep: one logical per registry variable, false when absent
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
% include: one entry per space, [] (all) when absent
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
% [ALPHA_SP,NB,BASE,BASIS_LIST,BSP,BIX] = BASIS_OPERATORS(S,INCL) lists the
% multi-indices of the included basis operators of every space.
% The basis operators Z_alpha of Sec. 9.1, per space: alpha_sp{k} holds     % MMP, 09/30/2026
% one multi-index per row over space k's own variables, nb(k) rows; base,   % MMP, 09/30/2026
% bsp, bix map global basis c to (space, row) and back; basis_list is the   % MMP, 09/30/2026
% output of the same name.                                                  % MMP, 09/30/2026
% INPUT  S:          'struct' of PARSE_INPUTS with field sep (1 x nv
%                    'logical'); read are M, nv, own, sep;
%        incl:       1 x M 'cell' of PARSE_OPTIONS.
% OUTPUT alpha_sp:   1 x M 'cell', alpha_sp{k} nb(k) x numel(own{k})
%                    'double';
%        nb:         1 x M 'double';
%        base:       1 x (M+1) 'double', [0, cumsum(nb)]: row i of space k
%                    is basis operator base(k)+i;
%        basis_list: N x (1+nv) 'double', N = sum(nb);
%        bsp, bix:   N x 1 'double', the space of basis operator c and its
%                    row in alpha_sp.
M = S.M;    nv = S.nv;  own = S.own;    sep = S.sep;                        % MMP, 09/30/2026

% % % Basis operators of each space. The multi-index runs over the space's
% OWN variables only: a variable the space does not have carries no
% indicator, so there is nothing to choose in that direction.
alpha_sp = cell(1,M);   nb = zeros(1,M);
for k = 1:M
    % Admissible entries per variable of the space: 1, 4 where sep is set,
    % else 1, 2, 3
    nk = numel(own{k});
    vals = cell(1,nk);
    for t = 1:nk
        if sep(own{k}(t)),  vals{t} = [1,4];  else,  vals{t} = [1,2,3];  end
    end
%   A = alpha_grid(vals);                                                   % MMP, 09/30/2026 (was)
    % Shared enumeration, direction 1 fastest (Sec. 9.1 multi-indices): a   % MMP, 09/30/2026
    % linear 'include' index means the same here and in 'sopquadvar'.       % MMP, 09/30/2026
    A = multiindex_grid(vals);                                              % MMP, 09/30/2026
    % The included ones: at least one, all distinct
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
% DEG_LIST = DEGREE_LIST(DEG,S,NB,BSP,BIX) expands the degree caps DEG of
% COPQUADVAR into one specification per basis operator.
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
% INPUT  deg:          the DEG of COPQUADVAR;
%        S:            'struct' of PARSE_INPUTS; read are M, nv, own;
%        nb, bsp, bix: 'double', of BASIS_OPERATORS.
% OUTPUT deg_list:     1 x N 'cell' of 'struct' with fields int, mult
%                      (1 x nv 'double'), joint ('double') and subset
%                      (1 x 2^(2*nv) 'double', or []).
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
    % The specification of basis operator c: that of its space, or its
    % entry in a per-basis cell
    k = bsp(c);     spec = deg_sp{k};
    if iscell(spec)
        if numel(spec)~=nb(k)
            error("A cell degree specification for space "+num2str(k)...
                  +" should have one entry per basis operator of that space.")
        end
        spec = spec{bix(c)};
    end
    % Expanded to int, mult (1 x nv), joint and subset
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
% [Z1C,Z2C,MDIM] = MONOMIAL_BASES(S,ALPHA_SP,DEG_LIST,BSP,BIX) returns the
% monomial vector Z^c of every basis operator c, within its degree caps:
% Z1c{c}, Z2c{c}: the monomials of basis c in (theta,s) and (theta,s'),     % MMP, 09/30/2026
% the two arguments 'sosquadvar' takes; mdim(c) its space dimension.        % MMP, 09/30/2026
% INPUT  S:        'struct' of PARSE_INPUTS; read are own, vars, vars_int,
%                  vars_dum, dims;
%        alpha_sp: 1 x M 'cell'; bsp, bix: N x 1 'double', of
%                  BASIS_OPERATORS;
%        deg_list: 1 x N 'cell' of DEGREE_LIST.
% OUTPUT Z1c, Z2c: 1 x N 'cell' of 'polynomial' column vectors, the
%                  monomials over [theta_1..theta_nv, s^k] and over
%                  [theta_1..theta_nv, s'^k];
%        mdim:     N x 1 'double', mdim(c) = m_k of the space of c.
own = S.own;    vars = S.vars;  vars_int = S.vars_int;                      % MMP, 09/30/2026
vars_dum = S.vars_dum;  dims = S.dims;  N = numel(bsp);                     % MMP, 09/30/2026

% Exponents are stored over [theta_1..theta_nv, s^k], so Z1c{c} and Z2c{c}
% differ only in the name of the second group.
Z1c = cell(1,N);    Z2c = cell(1,N);    mdim = zeros(N,1);
for c = 1:N
    % Caps in theta_1..theta_nv and in the variables of s^k
    k = bsp(c);     ok = own{k};
    caps_int  = deg_list{c}.int;
    caps_mult = deg_list{c}.mult(ok);
    % A multiplier in direction d identifies s_d with theta_d, so any s_d
    % dependence of the basis would be redundant.
    caps_mult(alpha_sp{k}(bix(c),:)==1) = 0;
    % Exponent grid; the monomial vectors in (theta,s) and (theta,s')
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gfun = psatz_weight(options,S)                                     % MMP, 09/30/2026
% GFUN = PSATZ_WEIGHT(OPTIONS,S) returns the weight g(theta) of
% OPTIONS.psatz:
% The weight g(theta) of options.psatz (0: g = 1), see the help.            % MMP, 09/30/2026
% INPUT  options: 'struct' with field psatz, a 'double': 0, 1, 2d+1 or
%                 2d+2;
%        S:       'struct' of PARSE_INPUTS; read are nv, vars_int, dom.
% OUTPUT gfun:    'polynomial' in the theta variables vars_int.
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
% [TERMS,NT] = BLOCK_TERMS(K,L,S,B) returns the terms of block (K,L) of
% Pop, from space L to space K:
% The terms Z_i* Q_ij Z_j (Sec. 9.2) of block (k,l), i over the basis       % MMP, 09/30/2026
% operators of space k and j over those of space l, integration variable    % MMP, 09/30/2026
% eliminated: terms{1:nt}, one sdopvar per substitution pair of each        % MMP, 09/30/2026
% (i,j), all on the one list B.Zd. HOT: runs the whole elimination. It      % MMP, 09/30/2026
% must not assign into Zd, Nprod or prog: writing one element of a          % MMP, 09/30/2026
% shared q-length cell copies all of it (27.5 ms at q = 2.7e6, measured).   % MMP, 09/30/2026
% INPUT  k, l:  'double', the output and the input space;
%        S:     'struct' of PARSE_INPUTS with field sep;
%        B:     'struct' with fields alpha_sp, nb, base (of
%               BASIS_OPERATORS), Nprod (N x N 'cell' of 'dpvar',
%               Nprod{c,e} = (I kron Z1c{c})' Q_ce (I kron Z2c{e})), Zd,
%               dmap, ndec (of COMMON_DVARS) and gfun (of PSATZ_WEIGHT).
% OUTPUT terms: 'cell', terms{1:nt} 'sdopvar' objects from space l to space
%               k, one per pair (i,j) and per substitution of 2 and 3 for
%               the entries 4 of the multi-indices in the shared
%               directions;
%        nt:    'double'.
nv = S.nv;  vars = S.vars;  vars_int = S.vars_int;  vars_dum = S.vars_dum;  % MMP, 09/30/2026
own = S.own;  dims = S.dims;  mask = S.mask;  dom = S.dom;  sep = S.sep;    % MMP, 09/30/2026
alpha_sp = B.alpha_sp;  nb = B.nb;  base = B.base;  Nprod = B.Nprod;        % MMP, 09/30/2026
gfun = B.gfun;  Zd = B.Zd;  dmap = B.dmap;  ndec = B.ndec;                  % MMP, 09/30/2026

% Adjoining a basis operator flips its lower and upper integrals; the
% full-domain integral is self-adjoint.
negmap = [1,3,2,4];

    % Variables and dimensions of the output space k and the input space l
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

    % Variable split (s,theta | s') of each product; the variables and
    % domains of the block
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

        % Monomial groups of the rows (s, theta) and the columns (s'), and
        % their sizes
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

        % Exponents of the shared directions laid out per variable of each
        % space
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

        % The merge maps over the matrix rows and columns
        Lmat = kron(speye(mk),ML).';
        Rmat = kron(speye(ml),MR);

        % % % Parameter cell of the term over the gamma cells of Sec. 6.1   % MMP, 09/30/2026
        params = gamma_params(Cgam,g1b,NL3,g2b,NR3,nloc,Lmat,Rmat, ...
                              locB,ndec,n3);                                % MMP, 09/30/2026

        % The term, on the one decision variable list Zd
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
% PARAMS = GAMMA_PARAMS(CGAM,G1B,NL3,G2B,NR3,NLOC,LMAT,RMAT,LOCB,NDEC,N3)
% returns the coefficients Lmat*C_q*Rmat of the gamma cells q of one term:
% The sdopvar parameter cell of one term: params.A{q}, params.B{q} (the     % MMP, 09/30/2026
% constant and decision parts of gamma cell q, B over the global list of    % MMP, 09/30/2026
% ndec variables), as a 3 x ... x 3 cell over the n3 shared variables.      % MMP, 09/30/2026
% Cgam{q}: triplets from 'int_semisep' on the block's own nloc sheets;      % MMP, 09/30/2026
% locB: their rows in the global list; Lmat, Rmat: the basis merge maps.    % MMP, 09/30/2026
% INPUT  Cgam:     'cell' of 3^n3 'struct' with fields i, j, v, m, n, the
%                  triplets of C_q, g1b*NL3 x (1+nloc)*g2b*NR3: column
%                  block 1 the constant part, block 1+r decision variable r
%                  of the term;
%        g1b, NL3, g2b, NR3: 'double', the row count g1b*NL3 and the column
%                  count g2b*NR3 of each block of C_q;
%        nloc:     'double', the number of decision variables of the term;
%        Lmat:     sparse (m_k*NLf) x (g1b*NL3), Rmat: sparse (g2b*NR3) x
%                  (m_l*NRf), the maps onto the merged bases of NLf and NRf
%                  monomials;
%        locB:     nloc x 1 'double', the indices of those variables in Zd;
%        ndec:     'double', numel(Zd);
%        n3:       'double', the number of shared directions.
% OUTPUT params:   'struct' with fields A and B, 'cell' arrays 3 x ... x 3
%                  over the n3 shared directions: A{q} the vectorized
%                  constant part, sparse (m_k*NLf*m_l*NRf) x 1, and B{q}
%                  the decision part, sparse ndec x (m_k*NLf*m_l*NRf).

        % Nothing between 'pack_sheets' and here mixes sheets (each stays a
        % contiguous block of columns), so the block's own nloc sheets carry
        % the elimination and the rows are scattered onto the global basis
        % once, here, from triplets; ndec sheets would pad every intermediate.
        params = struct();
        params.A = cell(numel(Cgam),1);
        params.B = cell(numel(Cgam),1);
%       KLR = kron(Rmat.',Lmat).';                                          % MMP, 09/26/2026 (was)
        % 'unpack_sheets' applies Lmat, Rmat and the scatter onto Zd as one
        % index remap of the triplets of Cgam{q}; no Kronecker product is
        % formed (the KLR of 09/26/2026 was superseded the same day).
        for q = 1:numel(Cgam)
%           [Aq,Bq] = unpack_sheets(Cgam{q},g1b*NL3,g2b*NR3,ndec);          % MMP, 09/22/2026 (was)
%           [Aq,Bq] = unpack_sheets(Cgam{q},g1b*NL3,g2b*NR3,nloc);          % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
%           [params.A{q},Bout] = lr_multiply(Lmat,Aq,Bq,Rmat);              % MMP, 09/22/2026 % MMP, 09/26/2026 (was)
%           [params.A{q},Bout] = lr_multiply(Lmat,Aq,Bq,Rmat,KLR);          % MMP, 09/26/2026 (was)
%           [params.A{q},params.B{q}] = lr_multiply(Lmat,Aq,Bq,Rmat);       % MMP, 09/22/2026 (was)
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
% Z = BY_SPACE(SRC,IDX,SEL) picks the entries IDX(SEL) of SRC:
% Exponent lists laid out per variable of one space: Z{t} = src{idx(t)}     % MMP, 09/30/2026
% where sel(t), and 0 (degree 0 only) elsewhere. Replaces four identical    % MMP, 09/30/2026
% scatter loops of 'block_terms' (initial coding 09/21/2026).               % MMP, 09/30/2026
% INPUT  src: 'cell' of exponent column vectors;
%        idx: 1 x n 'double', indices into src where sel is true;
%        sel: 1 x n 'logical'.
% OUTPUT Z:   1 x n 'cell'.

Z = repmat({0},1,numel(idx));                                               % MMP, 09/30/2026
Z(sel) = src(idx(sel));                                                     % MMP, 09/30/2026

end
% % % END change MMP, 09/30/2026 (split)                                    % MMP, 09/30/2026


% MMP, 09/30/2026: deleted the local function 'alpha_grid' (initial coding  % MMP, 09/30/2026
% 09/21/2026, no later markers), now 'multiindex_grid' (misc/conventions),  % MMP, 09/30/2026
% whose body is alpha_grid's.                                               % MMP, 09/30/2026


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = alpha_select(A,incl,nd,k)
% A = ALPHA_SELECT(A,INCL,ND,K) returns the rows of A that the include
% entry INCL of space K selects.
% INPUT  A:    n x nd 'double', the admissible multi-indices of space k,
%              one per row;
%        incl: [] (all rows); a 'logical' mask over the rows; an array with
%              nd columns, each row a multi-index among the rows of A; or
%              linear indices into the rows;
%        nd:   'double', the number of variables of space k;
%        k:    'double', the index of the space, for the error messages.
% OUTPUT A:    the selected rows; INCL itself for an array of
%              multi-indices.

% [] selects every row
if isempty(incl)
    return
end
% A logical mask, one entry per row
if islogical(incl)
    if numel(incl)~=size(A,1)
        error("A logical 'include' for space "+num2str(k)+" has one entry "...
              +"per basis operator of that space.")
    end
    A = A(incl(:),:);
% Rows of multi-indices. Validating them against the generated list rather
% than a numeric range is what makes a separable direction, which admits
% only 1 and 4, reject 2 and 3.
elseif nd>0 && size(incl,2)==nd
    if ~all(ismember(incl,A,'rows'))
        error("'include' for space "+num2str(k)+" contains an inadmissible "...
              +"multi-index.")
    end
    A = incl;
% Linear indices into the rows
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
% [C2,ZP,NL1,NR1] = INT_ONESIDED(C,G1,G2,NSHEET,ZTH,CLS,ONEIDX,DOM)
% integrates out the theta directions that are NOT shared by the block
% pair, the monomial theta_d^p in direction d becoming
%
%   cls(d)==1   int I_beta(s_d-theta_d) theta_d^p dtheta_d,  beta = oneidx(d)
%   cls(d)==2   int I_alpha(theta_d-s'_d) theta_d^p dtheta_d
%   cls(d)==0   int_a^b theta_d^p dtheta_d
%
% a polynomial in s_d, in s'_d, or a constant (alpha = oneidx(d) for
% cls(d)==2).
%
% INPUT  C:      (g1*NG) x ((1+nsheet)*g2) sparse, the packed coefficient
%                matrix, rows indexed (g1, theta) with theta fastest and
%                direction 1 slowest within it, columns (sheet, g2), NG the
%                number of theta monomials;
%        g1, g2: 'double', the row and column counts before the theta
%                monomials;
%        nsheet: 'double', the number of decision-variable sheets, sheet 1
%                (the constant part) not counted;
%        Zth:    1 x nv 'cell' of theta exponent column vectors;
%        cls:    1 x nv 'double', the direction classes: 3 shared (left
%                alone), 1 output space only, 2 input space only, 0
%                neither;
%        oneidx: 1 x nv 'double', multi-index entries, read only where cls
%                is 1 or 2;
%        dom:    nv x 2 'double', the domains.
% OUTPUT C2:     (g1*NL1*NG3) x ((1+nsheet)*g2*NR1) sparse, rows now
%                indexed (g1, new s monomials, remaining theta) and columns
%                (sheet, g2, new s' monomials), which is the layout
%                'int_semisep' consumes with g1*NL1 and g2*NR1 as its
%                dimensions, NG3 the number of remaining theta monomials;
%        Zp:     1 x nv 'cell', Zp{d} the new exponents produced in
%                direction d and [] where none are;
%        NL1, NR1: 'double', the sizes of the two new monomial groups.
%
% Sharing is what makes the semiseparable integral necessary: with a variable
% on both sides the result carries an indicator in s-s' and splits over
% gamma. With the variable on one side only, or on neither, the integral is
% an antiderivative of a monomial evaluated at one limit. Doing these first
% leaves 'int_semisep' only the shared directions.

% Per direction, the map W{d} from the theta exponents to the exponents of
% the result: the identity (shared), a definite integral (neither) or a
% one-sided integral
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

% Sizes: theta monomials; new s, new s' and remaining theta monomials
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

% Entries of T grouped by their theta monomial (column of T)
cnt = accumarray(tc,1,[NG,1]);
ptr = [0;cumsum(cnt)];
[~,ord] = sort(tc);
iLs = iL(ord);  iRs = iR(ord);  i3s = i3(ord);  tws = tw(ord);

% Fan-out of each nonzero of C; an all-zero C2 when there is no output
cc  = cnt(ith);
tot = sum(cc);
if tot==0
    C2 = sparse(g1*NL1*NG3,(1+nsheet)*g2*NR1);
    return
end
% Every 'repelem' result is forced to a column: 'repelem' of a 1 x 1 first
% argument returns a ROW, and 'pos' below is a column, so with exactly one
% nonzero in C the subtraction would expand into a tot x tot matrix, every
% entry emitted 'cc' times at the same (row,column) and summed by 'sparse':
% a block scaled by cc, with no error (nnz(C)==1 and a non-shared direction
% with fan-out above 1 only; invisible to a span comparison).
off = cumsum([0;cc(1:end-1)]);
pos = (1:tot).' - colv(repelem(off,cc));        % 1..cc(u) within nonzero u
sel = colv(repelem(ptr(ith),cc)) + pos;

% Row, column and value of every output entry
newrow = (colv(repelem(arow,cc))*NL1 + iLs(sel))*NG3 + i3s(sel) + 1;
newcol = colv(repelem(sh,cc))*(g2*NR1) + colv(repelem(jc,cc))*NR1 + iRs(sel) + 1;
newval = colv(repelem(v,cc)).*tws(sel);

C2 = sparse(newrow,newcol,newval,g1*NL1*NG3,(1+nsheet)*g2*NR1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = colv(x)
% X = COLV(X) returns X as a column, so that an arithmetic combination of
% two such vectors cannot silently become a matrix by implicit expansion.
% INPUT  x: array.
% OUTPUT x: numel(x) x 1, the entries of X in linear order.

x = reshape(x,[],1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Zo,W] = onesided_int(P,side,idx,lim)
% [ZO,W] = ONESIDED_INT(P,SIDE,IDX,LIM) integrates theta^p against a
% single indicator over the domain, returning the exponents Zo of the
% result in the surviving spatial variable and the map W from the theta
% exponents P to them.
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
%
% INPUT  P:    'double' column of the nonnegative exponents p of theta_d;
%        side: 'double', 1 or 2 as above;
%        idx:  'double', the multi-index entry: 1 multiplier (theta_d
%              becomes the variable), 2 or 3 as above, 4 full-domain
%              integral;
%        lim:  1 x 2 'double', the domain [a,b].
% OUTPUT Zo:   'double' column of the exponents of the variable (0 alone
%              for idx 4);
%        W:    numel(Zo) x numel(P) sparse, column t the coefficients of
%              the integral of theta^P(t) on the monomials of Zo.

% Multiplier: substitution; full-domain integral: definite integral
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

% Lower or upper integral: whether the variable is the upper limit
a = lim(1);     b = lim(2);     P = reshape(P,[],1);    np = numel(P);
var_is_upper = (side==1 && idx==2) || (side==2 && idx==3);

% Exponents of the result; the constant and variable coefficients of each
% antiderivative, in the rows of Zo of exponent 0 and p+1
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
% [ZO,W] = DEFINITE_INT(P,LIM) returns int_a^b theta^p dtheta for each
% exponent p in P: a single constant row.
% INPUT  P:   'double' vector of the nonnegative exponents;
%        lim: 1 x 2 'double', the domain [a,b].
% OUTPUT Zo:  0, the exponent of the constant;
%        W:   1 x numel(P) sparse, W(t) = (b^(p+1)-a^(p+1))/(p+1) for
%             p = P(t).

% One coefficient per exponent, on the constant monomial
a = lim(1);     b = lim(2);     P = reshape(P,[],1);    np = numel(P);
Zo = 0;
w = zeros(1,np);
for t = 1:np
    p = P(t);
    w(t) = (b^(p+1)-a^(p+1))/(p+1);
end
W = sparse(ones(1,np),1:np,w,1,np);

end
