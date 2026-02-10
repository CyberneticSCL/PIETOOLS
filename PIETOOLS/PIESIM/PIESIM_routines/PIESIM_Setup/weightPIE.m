%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% weightPIE.m     PIETOOLS 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PIE = weightPIE(PIE,k)
% Multiplies the left-hand side of the PIE equation by s^k and modifies the
% corresponding T, Tu, Tw and mapping operators. Utilized in cylindrical
% coordinates examples.
% Inputs:
% 1) PIE - a PIE structure
% 2) k - weight
% Both inputs are REQUIRED
% Outputs: PIE - a PIE with the left-hand side weighted by s^k, such that
% s^k T dot(x_f) + s^k Tw dot(w) + s^k Tu dot(u) = A x_f + B1 w + B2 u
% x = Tumap u + Twmap w + Tmap x
% Operators Tumap, Twmap, and Tmap are stored in the PIE.misc construct,
% while the weighted T, Tw and Tu operators are stored in the standard
% PIE constuct

% Initial coding: YP 1/7/2026

<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/PIESIM_Setup/weightPIE.m
pvar s1;

=======
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/PIESIM_Setup/weightPIE.m
        % Initialize mapping operators to be unmodified PIE operators

        PIE.misc.Tmap = PIE.T;
        PIE.misc.Tumap = PIE.Tu;
        PIE.misc.Twmap = PIE.Tw;

        % Modify PIE operators that describe dynamics due to the multiplication
        % of th PDE by r^k (k=weight)
    
        opvar W;
        [m, n] = size(PIE.T.R.R0);
        W.dim = [0 0; m n];
        W.I = PIE.T.I;
        W.var1 = s1;
        W.R.R0 = s1^k * eye(m);
        W.R.R1 = 0;
        W.R.R2 = 0;
        
        PIE.T = W * PIE.T;
        PIE.Tu = W * PIE.Tu;
        PIE.Tw = W* PIE.Tw;
<<<<<<< Updated upstream:PIETOOLS/PIESIM/PIESIM_routines/PIESIM_Setup/weightPIE.m
end
=======
end
>>>>>>> Stashed changes:PIESIM/PIESIM_routines/PIESIM_Setup/weightPIE.m
