function [prog,Pop,Qcell,basis_list] = poscopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,POP,QCELL,BASIS_LIST] = POSCOPVAR(PROG,DIMS,SPACES,DOM,DEG,OPTIONS)
% declares in the LPI program PROG a positive semidefinite, self-adjoint
% 'cdopvar' decision operator on a concatenation of spaces, parameterized
% as in Sec. 8.4 of the sopvar document,
%
%       Pop: X -> X,     X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M],
%
%       Pop = sum_{b,b'} Z_b'*(g*Q_{b,b'})*Z_{b'},     Q >= 0,
%
% with ONE Gram matrix Q of decision variables over every pair (b,b') of
% basis operators of all spaces, so that Pop is positive as a whole, and
% g(theta) >= 0 the weight of OPTIONS.psatz. An empty s^k is the space
% R^{m_k}. Basis operator b of space k carries a multi-index alpha over the
% variables of s^k and maps L_2^{m_k}[s^k] to L_2[theta], theta over every
% variable of the registry:
%
%       (Z_b x)(theta) = int I_alpha(theta-s) ...
%                            (I_{m_k} kron Z^alpha(theta,s)) x(s) ds,
%
% I_alpha(r) = prod_d I_{alpha_d}(r_d) over the variables of s^k, alpha_d =
% 1 the multiplier (delta), 2 the lower integral (s_d <= theta_d), 3 the
% upper integral (s_d >= theta_d), 4 the full integral; Z^alpha(theta,s)
% is the vector of monomials of the degrees DEG, and on R^{m_k},
% Z_b x = (I_{m_k} kron Z^alpha(theta)) x. Pop is also appended to the
% cell prog.sop_pos, from which LPI_INEQ_SOP sizes an inequality that
% depends on Pop.
%
% INPUT
% - prog:        'struct', the LPI program (LPIPROGRAM_SOP);
% - dims:        double, M x 1 positive integers, dims(k) = m_k the number
%                of components of space k; a scalar is used for every
%                space; empty for no space (POP is then 0 x 0, with no
%                decision variable, and PROG is returned unchanged);
% - spaces:      1 x M cell of 'cellstr', spaces{k} the names of the
%                variables of s^k ({} for R^{m_k}); a plain 'cellstr' is
%                one space, so several spaces are nested, as {{},{'s'}};
% - dom:         double, nv x 2, dom(d,:) the interval of variable d of
%                the registry (the sorted union of all spaces), or one
%                1 x 2 row for every variable; or a 'struct' with fields
%                vars ('cellstr') and dom (double, one row per name);
% - deg:         the degrees of Z^alpha(theta,s): a scalar (int and mult),
%                or a 'struct' with the optional fields
%                  int     scalar or 1 x nv, the degree in each theta_d
%                          (default 1);
%                  mult    scalar or 1 x nv, the degree in each s_d of
%                          s^k, 0 where alpha_d = 1 (default int);
%                  joint   scalar, a cap on the total degree (default
%                          none);
%                  subset  2^(2nv) caps, one per subset of [theta, s^k],
%                          for a space with every registry variable only;
%                or a 1 x M cell of these, one per space, each of which
%                may be a cell with one per basis operator of the space;
% - options:     (optional) 'struct' with the optional fields
%                  type     char, 'pos' (any other value is an error);
%                  psatz    double, the weight g: 0 (default, g = 1), 1
%                           (g = prod_d (theta_d-a_d)*(b_d-theta_d)), or
%                           2d+1 / 2d+2, d = 1..nv (g = (theta_d-a_d)/L_d
%                           or (b_d-theta_d)/L_d, L_d = b_d-a_d); any
%                           other value is an error;
%                  sep      logical, scalar or 1 x nv: where true, the
%                           lower and upper basis operators of that
%                           direction are replaced by the full integral
%                           (alpha_d = 4);
%                  include  1 x M cell, entry k selecting basis operators
%                           of space k: multi-index rows over the
%                           variables of s^k, a logical mask, or linear
%                           indices (default all);
%                  path     char, 'fast' (default) or 'general', the
%                           construction of the map from Q to the kernels
%                           of Pop (the same Pop up to rounding);
%                  cache    logical (default true): keep the position
%                           maps of the session; POSCOPVAR_DIRECT('clear')
%                           empties them.
%
% OUTPUT
% - prog:        'struct', the program with Q declared positive
%                semidefinite and Pop appended to prog.sop_pos;
% - Pop:         'cdopvar' object, M x M blocks on X, the operator above;
% - Qcell:       N_b x N_b cell, Qcell{b,b'} the 'cellstr' names of the
%                decision variables of block (b,b') of Q, N_b the number
%                of basis operators of all spaces;
% - basis_list:  double, N_b x (1+nv), row b the basis operator b: its
%                space k, then its multi-index entry alpha_d in each
%                registry direction, 0 where s^k lacks variable d; rows
%                ordered by space, as Qcell.
%
% See also POSCOPVAR_DIRECT, COPQUADVAR, POSSOPVAR, SOPQUADVAR, POSLPIVAR,
%          LPI_EQ_CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - poscopvar
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
%                  mopquadvar -> copquadvar, posmopvar -> poscopvar. File was
%                  'posmopvar.m'.
% MMP, 10/17/2026: Built by 'poscopvar_direct' instead of 'copquadvar': the
%                  same operator and decision variables (24 cases to
%                  4.4e-16, test_poscopvar_direct) at 2.0-3.8x less time
%                  cold and 7-19x warm, and in 3-D a third of the peak
%                  memory (17 operators, 10/16/2026); in the executives the
%                  storage builds 5-15x faster, the same program and values
%                  to 10 digits. Cost: the position map of each distinct
%                  shape stays cached for the session, about 16 B per
%                  coefficient nonzero (117 MB for 3-D at int = mult = 1,
%                  q = 5.0e5); poscopvar_direct('clear') releases it,
%                  options.cache = false skips it. 'copquadvar' stays the
%                  routine for type 'sym' and the reference of the tests.
%                  The case with every space R^q (no spatial variable, the
%                  W of the H2 Q-form executives) needs poscopvar_direct's
%                  fix of the same date; before it, pos_map errored there
%                  (found by the chain test on PIETOOLS_H2_norm_c_sop).
%                  options.psatz keeps copquadvar's check (one code, 0, 1
%                  or a face), since poscopvar_direct accepts more (a row
%                  of codes, the box quadratics 2nv+2+d) and
%                  test_copquadvar_faces asserts the rejections.
% MMP, 10/21/2026: Pop is recorded in prog.sop_pos, from which lpi_ineq_sop
%                  sizes the certificate of an inequality that depends on
%                  Pop (maintainer: lpi_ineq_sop needs no hand-set options).
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                  inputs, outputs), step comments in the body; code
%                  unchanged. The deleted header text was the wrapper
%                  description (its earlier version kept as (was) lines),
%                  the references to copquadvar and poscopvar_direct for
%                  the inputs and outputs, a test name and the reason
%                  options.type is refused. The help text of the
%                  10/21/2026 entry is in the paragraph and the OUTPUT
%                  entry of prog, without markers.
% MMP, 10/28/2026: Empty DIMS (no space: an empty channel of a PIE, as the
%                  W of the H2 executives on R^0 when the PIE has no
%                  disturbance or no regulated output) returns the 0 x 0
%                  operator with no decision variable. Before, the space
%                  parser read the empty list as one R^q space and failed
%                  (maintainer: empty channels should not raise an error).

% Check the inputs; OPTIONS a struct, of type 'pos' only
if nargin<5
    error("Not enough input arguments.")
end
if nargin<6 || isempty(options)
    options = struct();
end
if ~isa(options,'struct')
    error("Options should be specified as a 'struct' object.")
end
if isfield(options,'type') && ~isempty(options.type) ...
        && ~strcmp(char(options.type),'pos')
    error("'poscopvar' declares a positive operator; for type '%s' call "...
          +"'copquadvar' directly.",char(options.type))
end
options.type = 'pos';

% No space (an empty channel of a PIE, as the W of an H2 LPI on R^0): the   % MMP, 10/28/2026
% 0 x 0 operator, with no block and no decision variable; PROG unchanged    % MMP, 10/28/2026
if isempty(dims)                                                            % MMP, 10/28/2026
    meta0 = struct('vars',{cell(1,0)},'dom',zeros(0,2),'space_out',false(0,0),...
                   'space_in',false(0,0),'dim_out',zeros(0,1),'dim_in',zeros(0,1),...
                   'Zd',{cell(0,1)});                                       % MMP, 10/28/2026
    Pop = cdopvar(cell(0,0),meta0);     Qcell = {};     basis_list = zeros(0,1); % MMP, 10/28/2026
    return                                                                  % MMP, 10/28/2026
end                                                                         % MMP, 10/28/2026

% [prog,Pop,Qcell,basis_list] = copquadvar(prog,dims,spaces,dom,deg,options); % MMP, 10/17/2026 (was)
meta = parse_copvar_spaces(dims,spaces,dom);    % the registry, for nv      % MMP, 10/17/2026
% copquadvar's contract for options.psatz is kept: one code, 0, 1 or a face % MMP, 10/17/2026
% 2d+1 / 2d+2. poscopvar_direct also takes a row of codes and the box       % MMP, 10/17/2026
% quadratics 2nv+2+d (the terms of lpi_ineq_sop); poscopvar rejects them,   % MMP, 10/17/2026
% as before (test_copquadvar_faces part (i)), with copquadvar's message.    % MMP, 10/17/2026
nv = numel(meta.vars);      ps = 0;                                         % MMP, 10/17/2026
if isfield(options,'psatz') && ~isempty(options.psatz)                      % MMP, 10/17/2026
    ps = options.psatz;                                                     % MMP, 10/17/2026
end                                                                         % MMP, 10/17/2026
if ~isscalar(ps) || ~ismember(ps,[0,1,3:2*nv+2])                            % MMP, 10/17/2026
    error("'psatz' should be 0, 1, or a face code 2d+1 / 2d+2 with d "...
          +"in 1..nv, here nv = "+num2str(nv)+"; 2 is "...
          +"poslpivar_2d's ball, not the box.")                             % MMP, 10/17/2026
end                                                                         % MMP, 10/17/2026
% Declare Q and build Pop
[prog,Pop,Qcell,info] = poscopvar_direct(prog,dims,spaces,dom,deg,options); % MMP, 10/17/2026
basis_list = info.basis_list;                   % as copquadvar's           % MMP, 10/17/2026
% The record of lpi_ineq_sop: an inequality that depends on Pop is sized    % MMP, 10/21/2026
% from Pop's support, not from every value of its decision variables.       % MMP, 10/21/2026
% MATLAB copies on write, so the record shares Pop's arrays.                % MMP, 10/21/2026
if ~isfield(prog,'sop_pos'),    prog.sop_pos = {};  end                     % MMP, 10/21/2026
prog.sop_pos{end+1} = Pop;                                                  % MMP, 10/21/2026

end
