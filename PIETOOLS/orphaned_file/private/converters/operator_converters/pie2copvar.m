function PIEc = pie2copvar(PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIEC = PIE2COPVAR(PIE) takes a 'pie_struct' whose operators are 'opvar'
% (1-D) or 'opvar2d' (2-D) objects and returns the same 'pie_struct' with
% each operator a 'copvar' container of 'sopvar' blocks.
%
% INPUTS
% - PIE:    'pie_struct' of the PIE
%
%               T*x_t + Tw*w_t + Tu*u_t = A*x  + B1*w  + B2*u,
%                                     z = C1*x + D11*w + D12*u,
%                                     y = C2*x + D21*w + D22*u,
%
%           the twelve operators 'opvar' or 'opvar2d' objects; absent
%           operators become zero operators of the dimensions the other
%           operators imply;
%
% OUTPUTS
% - PIEc:   the initialized 'pie_struct', with the vars, dom, component
%           tables and misc of PIE, and each operator (K included, if
%           present) a 'copvar' container whose block (i,j) maps the j-th
%           input space of nonzero dimension to the i-th output space of
%           nonzero dimension, the spaces in the order (R, L2[s]) in 1-D
%           and (R, L2[s1], L2[s2], L2[s1,s2]) in 2-D.
%
% NOTES
% A zero operator between spaces of nonzero dimension (D11 = 0, Tw = 0)
% converts to zero blocks, which keep its spaces defined. An operator with
% no input or no output space (an empty channel) converts to a 'copvar'
% without spaces on that side.
%
% A PIE whose T is already a 'copvar' is returned as is.
%
% See also OP2COPVAR_SOP, OPVAR2COPVAR, COPVAR2OPVAR, COPVAR2OPVAR2D,
% PIE_STRUCT.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie2copvar
%
% Copyright (C)2026  PIETOOLS Team
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
% Initial coding MMP, 10/27/2026

% Check the input
if ~isa(PIE,'pie_struct')
    error('pie2copvar:pie','The PIE should be a ''pie_struct''.');
end
PIEc = PIE;
if isa(PIE.T,'copvar'),     return,     end         % already converted
% Initialize: absent operators become zero operators of the dimensions the
% other operators imply; dimensions, domains and variables are checked; the
% operators are converted to copvar objects, an operator with an empty
% channel (no input or no output space) to a copvar object without spaces
% on that side.
PIEc = initialize(PIE);
end
