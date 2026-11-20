%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weightPIE.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PIE = weightPIE(PIE,k)
% Multiplies the left-hand side of the PIE equation by s^k and modifies the
% corresponding T, Tu, Tw and mapping operators.
% Inputs:
% 1) PIE - a PIE structure, whose operators T, Tu, Tw are 'struct' objects
%    of kernels (fields dim, Q2, R.R0, R.R1, R.R2) in the spatial variable
%    s = s1, with no ODE state (T.dim(1,1) = 0, else an error is thrown)
% 2) k - weight
% Both inputs are REQUIRED
% Outputs: PIE - a PIE with the left-hand side weighted by s^k, such that
% s^k T dot(x_f) + s^k Tw dot(w) + s^k Tu dot(u) = A x_f + B1 w + B2 u
% x = Tumap u + Twmap w + Tmap x
% Operators Tumap, Twmap, and Tmap are stored in the PIE.misc construct,
% while the weighted T, Tw and Tu operators are stored in the standard
% PIE constuct
% (weighting multiplies the kernels Q2, R0, R1, R2 of each operator by s^k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - weightPIE
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding: YP 1/7/2026

% Spatial variable s1 of the weight
pvar s1;

        % Initialize mapping operators to be unmodified PIE operators

        pvar s1;

        PIE.misc.Tmap = PIE.T;
        PIE.misc.Tumap = PIE.Tu;
        PIE.misc.Twmap = PIE.Tw;

        % Modify PIE operators that describe dynamics due to the multiplication
        % of the PDE by r^k (k=weight)
    
        % W*X for W = s1^k I on L2: Q2, R0, R1, R2 of X times s1^k; X must
        % have no R^n output, since W has no R^n input
        for op = {'T','Tu','Tw'}
            X = PIE.(op{1});
            if X.dim(1,1)~=0
                error("The weight requires a PIE without ODE states.");
            end
            X.Q2 = s1^k*X.Q2;
            X.R.R0 = s1^k*X.R.R0;
            X.R.R1 = s1^k*X.R.R1;
            X.R.R2 = s1^k*X.R.R2;
            PIE.(op{1}) = X;
        end
end
