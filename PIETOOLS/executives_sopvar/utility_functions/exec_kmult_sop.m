function [prog,X] = exec_kmult_sop(prog,X,K,T,dom,like)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PROG,X] = EXEC_KMULT_SOP(PROG,X,K,T,DOM,LIKE) returns X + K'*Z + Z'*K,
% with Z a free operator added to PROG, from the state space of the PIE
% (the input space of T) to the output space of the constraint K.
%
% For a PIE whose state satisfies K x = 0, x'(K'*Z + Z'*K)x = 0 on ker K,
% so a negativity condition X + K'*Z + Z'*K <= 0 imposes X <= 0 on ker K
% only, as the dynamics require.
%
% INPUT
% - prog:   'struct', the LPI program;
% - X:      'cdopvar' object, the block of a negativity condition that acts
%           on the state, from and to the input space of T;
% - K:      'copvar' object, the constraint K x = 0 on the state;
% - T:      'copvar' object, the operator T of the PIE;
% - dom:    nv x 2 double, the domain of the spatial variables;
% - like:   'cdopvar' object, the storage operator whose degrees Z takes.
%
% OUTPUT
% - prog:   'struct', the program with the decision variables of Z;
% - X:      'cdopvar' object, X + K'*Z + Z'*K.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - exec_kmult_sop
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
% Initial coding MMP, 11/07/2026: shared by the primal LPIs of
%                PIETOOLS_Hinf_gain_sop, PIETOOLS_H2_norm_sop,
%                PIETOOLS_Hinf_estimator_sop and PIETOOLS_H2_estimator_sop
%                for PIEs with a constraint K (convert_PIETOOLS_PDE, e.g.
%                periodic boundary conditions).

% Z free, from the state space (input space of T) to the output space of K,
% with the degrees of the storage
[spK,dmK] = copvar_space_list(K,'out');     [spI,dmI] = copvar_space_list(T,'in');
[prog,Z] = lpivar_cdopvar(prog,struct('out',dmK,'in',dmI),struct('out',{spK},'in',{spI}),dom,get_lpivar_degs_sop(like));
% The multiplier, zero on ker K
X = X + K'*Z + Z'*K;
end
