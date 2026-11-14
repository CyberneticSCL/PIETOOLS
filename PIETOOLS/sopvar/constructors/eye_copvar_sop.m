function Iop = eye_copvar_sop(dims,spaces,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IOP = EYE_COPVAR_SOP(DIMS,SPACES,DOM) returns the identity operator on
%
%   L_2^{q_1}[s^1] x ... x L_2^{q_M}[s^M]
%
% as an M x M 'copvar': block (i,i) the multiplier I_{q_i} on space i, every
% other block zero. It is the container counterpart of
% mat2opvar(eye(n),dim,vars,dom), with any number of spatial variables.
%
% IOP = EYE_COPVAR_SOP(X) returns the identity on the output spaces of the
% 'copvar' or 'cdopvar' object X, with the domain of X: for an operator P
% of POSCOPVAR, P + e*EYE_COPVAR_SOP(P) is P + e*I.
%
% INPUT
% - dims:   M x 1 component counts q_i (a scalar is expanded), or a struct
%           with field 'out' (a field 'in' must then equal it);
% - spaces: 1 x M cell of cellstr variable names, {} being R^q, or a
%           struct with field 'out' (a field 'in' must then equal it);
% - dom:    as for 'mat2copvar_sop';
% - X:      'copvar' or 'cdopvar' object, in place of DIMS, SPACES and DOM
%           (second form above).
%
% OUTPUT
% - Iop:    M x M 'copvar', diagonal, passing 'verify'. A space listed twice
%           gives two separate components, and Iop is still the identity.
%
% NOTES
% Cost: O(M^2) block visits, O(q_i^2) per diagonal block, no decision
% variables.
%
% See also MAT2COPVAR_SOP, ZEROS_COPVAR_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eye_copvar_sop
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
% Initial coding MMP, 09/29/2026. Tier 1c: the Iw, Iz and eppos*I of the
%                executives, without opvar2copvar(mat2opvar(...)).
% MMP, 09/30/2026: Read the spaces with the shared 'parse_copvar_spaces'
%                (sopvar/misc/conventions); 'private/spaces2meta_sop', a
%                verbatim copy of the same parser, is deleted. Same meta on
%                every input; a non-cell space list now raises the parser's
%                message instead of MATLAB's brace-indexing error.
% MMP, 11/14/2026: Operator form EYE_COPVAR_SOP(X): the output spaces and
%                the domain of X, read by COPVAR_SPACE_LIST (maintainer,
%                11/14/2026: P + 1e-4*I restated the spaces of P and was too
%                complicated).

% Operator form EYE_COPVAR_SOP(X): the output spaces and domain of X        % MMP, 11/14/2026
if isa(dims,'copvar') || isa(dims,'cdopvar')                                % MMP, 11/14/2026
    narginchk(1,1);                                                         % MMP, 11/14/2026
    X = dims;   [spaces,dims] = copvar_space_list(X,'out');                 % MMP, 11/14/2026
    dom = struct('vars',{X.vars},'dom',X.dom);                              % MMP, 11/14/2026
end                                                                         % MMP, 11/14/2026
% meta = spaces2meta_sop(dims,spaces,dom);                                  % MMP, 09/30/2026 (was)
meta = parse_copvar_spaces(dims,spaces,dom);                                % MMP, 09/30/2026
if ~isequal(meta.space_out,meta.space_in) || ~isequal(meta.dim_out,meta.dim_in)
    error('eye_copvar_sop:notSquare',...
        'The identity needs the same spaces and dimensions in and out.')
end
Iop = mat2copvar_grid(speye(sum(meta.dim_out)),meta,true);

end
