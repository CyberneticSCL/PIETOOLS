function [prog,Pop,Qcell,basis_list] = posmopvar(prog,dims,spaces,dom,deg,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [prog,Pop,Qcell,basis_list] = POSMOPVAR(prog,dims,spaces,dom,deg,options)
% declares a positive semidefinite, self-adjoint 'mdopvar' decision operator
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
% This is a wrapper. The construction lives in 'mopquadvar', which assembles
% that quadratic form and takes the variable type as an argument;
% 'posmopvar' fixes the type to 'pos'. Every input and output, and every
% option other than 'type', is exactly as documented in 'mopquadvar' -- see
% there for dims, spaces, dom, deg, options.psatz, options.sep,
% options.include and the returned Qcell and basis_list.
%
% Passing options.type here is an error rather than a silent override, since
% a 'sym' variable would not be a positive operator and the name would then
% be wrong.
%
% See also MOPQUADVAR, POSSOPVAR, SOPQUADVAR, POSLPIVAR, LPI_EQ_MDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - posmopvar
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
    error("'posmopvar' declares a positive operator; for type '%s' call "...
          +"'mopquadvar' directly.",char(options.type))
end
options.type = 'pos';

[prog,Pop,Qcell,basis_list] = mopquadvar(prog,dims,spaces,dom,deg,options);

end
