function PIE = exec_pie_sop(PIE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIE = EXEC_PIE_SOP(PIE) prepares the PIE that a container executive
% (executives_sopvar) receives. A PIE whose T is not yet a 'copvar' object  % MMP, 11/09/2026
% is initialized (INITIALIZE). A container PIE, whose T is a                % MMP, 11/09/2026
% 'copvar' object (as PIE2COPVAR returns, or built in any number of
% spatial variables), is returned with every absent operator ([]) set to
% the double 0, so that the zero tests of the executives (Op==0, a logical
% scalar for a 'copvar' and for 0) hold for it as for an initialized PIE.
%
% INPUT
% - PIE:  'pie_struct' object of the PIE
%             T*x_t + Tw*w_t + Tu*u_t = A*x  + B1*w  + B2*u,
%                                   z = C1*x + D11*w + D12*u,
%                                   y = C2*x + D21*w + D22*u,
%         with T a 'copvar' object and each other operator a 'copvar'       % MMP, 11/09/2026
%         object or [] (absent), or not yet initialized; for a              % MMP, 11/09/2026
%         container PIE, vars (nv x 2 'polynomial', the spatial variables
%         and their dummies) and dom (nv x 2 double) define the dimension
%         nv = PIE.dim.
%
% OUTPUT
% - PIE:  'pie_struct' object: INITIALIZE(PIE) when T is not a 'copvar'; for a % MMP, 11/09/2026
%         container PIE, the input with each of T, Tw, Tu, A, B1, B2, C1,
%         D11, D12, C2, D21, D22 that is [] replaced by 0.                  % MMP, 11/09/2026
%         Every operator of the output is a 'copvar' object or 0.
%
% NOTES
% A PIE other than a 'pie_struct', and a container PIE whose vars and dom
% do not define a dimension nv >= 1, are errors.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - exec_pie_sop
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
% Initial coding MMP, 10/29/2026: the container executives take container
%                PIEs, in any number of spatial variables (maintainer: the
%                executives must support n variables, at least 3-D).
%                INITIALIZE errors on a container PIE, and an absent
%                operator [] makes Op==0 empty, which '&&' rejects.
% MMP, 11/04/2026: An 'opvar' or 'opvar2d' operator in a container PIE is
%                converted here. The executives no longer call
%                OP2COPVAR_SOP on the PIE operators, which converted such
%                an operator; this keeps that case working and makes every
%                operator of the output a 'copvar' or 0.
% MMP, 11/09/2026: Dead opvar branch removed (maintainer, 11/09/2026): the
%                conversion of an 'opvar'/'opvar2d' operator beside a
%                'copvar' T (11/04/2026) is gone with the retired classes;
%                help to match. INITIALIZE still serves a PIE whose T is
%                not a 'copvar'.

% A pie_struct only, as every executive requires
if ~isa(PIE,'pie_struct')
    error('exec_pie_sop:pie','The PIE should be a ''pie_struct'' object.');
end
% A PIE not yet initialized (T not a 'copvar' object): INITIALIZE           % MMP, 11/09/2026
if ~isa(PIE.T,'copvar')
    PIE = initialize(PIE);
    return
end
% A container PIE: its dimension from vars and dom (NaN when they differ)
if ~(PIE.dim>=1)
    error('exec_pie_sop:dim',['A container PIE needs vars (nv x 2 ''polynomial'') and dom ' ...
          '(nv x 2) of the same length nv >= 1.']);
end
% Absent operators to 0. Only [] is replaced: a 'copvar' with no space on a
% side (an empty channel) stays as it is.
ops = {'T','Tw','Tu','A','B1','B2','C1','D11','D12','C2','D21','D22'};
for k = 1:numel(ops)
    X = PIE.(ops{k});
    if isnumeric(X) && isempty(X),  PIE.(ops{k}) = 0;   end
%   if isa(X,'opvar') || isa(X,'opvar2d'),  PIE.(ops{k}) = op2copvar_sop(X);    end % MMP, 11/09/2026 (was); MMP, 11/04/2026
end
end
