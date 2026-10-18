function [prog,Pop,Qcell,basis_list] = poscopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,Qcell,basis_list] = POSCOPVAR(prog,dims,spaces,dom,deg,options)
% declares a positive semidefinite, self-adjoint 'cdopvar' decision operator
% on a concatenation of mixed spaces,
%
%       Pop: X -> X,     X = L_2^{m_1}[s^1] x ... x L_2^{m_M}[s^M],
%       Pop = Pop* >= 0,
%
% parameterized as in Sec. 8.4 of the sopvar document,
%
%       Pop = sum_{k,l} sum_{i,j} (Z_{alpha_i^k})* Q_{(ki),(lj)} (Z_{alpha_j^l})
%
% with a SINGLE Q>=0 spanning every (space, multi-index) pair, so that the
% container is positive as a whole and not merely blockwise. An empty space
% is R^{m_k}, which is what makes this cover the mixed finite/infinite case
% of a PIE that 'possopvar' does not.
%
% This is a wrapper. The construction lives in 'copquadvar', which assembles % MMP, 10/17/2026 (was)
% that quadratic form and takes the variable type as an argument;           % MMP, 10/17/2026 (was)
% 'poscopvar' fixes the type to 'pos'. Every input and output, and every    % MMP, 10/17/2026 (was)
% option other than 'type', is exactly as documented in 'copquadvar' -- see % MMP, 10/17/2026 (was)
% there for dims, spaces, dom, deg, options.psatz, options.sep,             % MMP, 10/17/2026 (was)
% options.include and the returned Qcell and basis_list.                    % MMP, 10/17/2026 (was)
% This is a wrapper of 'poscopvar_direct', which maps the Gram straight to  % MMP, 10/17/2026
% the kernel coefficients: the operator, decision variables (names and      % MMP, 10/17/2026
% order), Qcell and basis_list are those of 'copquadvar' with type 'pos'    % MMP, 10/17/2026
% (test_poscopvar_direct). Every input and output, and every option other   % MMP, 10/17/2026
% than 'type', is as documented in 'copquadvar' -- see there for dims,      % MMP, 10/17/2026
% spaces, dom, deg, options.psatz, options.sep, options.include and the     % MMP, 10/17/2026
% returned Qcell and basis_list; options.path and options.cache are those   % MMP, 10/17/2026
% of 'poscopvar_direct'. 'copquadvar' builds the same operator through the  % MMP, 10/17/2026
% pair integrals of its basis operators and takes type 'sym' as well.       % MMP, 10/17/2026
%
% Passing options.type here is an error rather than a silent override, since
% a 'sym' variable would not be a positive operator and the name would then
% be wrong.
%
% See also COPQUADVAR, POSSOPVAR, SOPQUADVAR, POSLPIVAR, LPI_EQ_CDOPVAR.    % MMP, 10/17/2026 (was)
% See also POSCOPVAR_DIRECT, COPQUADVAR, POSSOPVAR, SOPQUADVAR, POSLPIVAR,  % MMP, 10/17/2026
%          LPI_EQ_CDOPVAR.                                                  % MMP, 10/17/2026
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
[prog,Pop,Qcell,info] = poscopvar_direct(prog,dims,spaces,dom,deg,options); % MMP, 10/17/2026
basis_list = info.basis_list;                   % as copquadvar's           % MMP, 10/17/2026

end
