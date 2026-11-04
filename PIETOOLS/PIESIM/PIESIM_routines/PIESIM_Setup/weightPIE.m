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
% MMP, 11/04/2026: Weight the kernels directly instead of composing with an
%   opvar W: PIESIM now passes each operator as a struct of its kernels
%   (pie2piesim), which has no operator product. W = s1^k I on L2 is a
%   multiplier without R^n part, so W*X multiplies Q2, R0, R1 and R2 of X
%   by s1^k and needs X without R^n output, as the product did.

pvar s1;

        % Initialize mapping operators to be unmodified PIE operators

        pvar s1;

        PIE.misc.Tmap = PIE.T;
        PIE.misc.Tumap = PIE.Tu;
        PIE.misc.Twmap = PIE.Tw;

        % Modify PIE operators that describe dynamics due to the multiplication
        % of the PDE by r^k (k=weight)
    
%       opvar W;                                                            % MMP, 11/04/2026 (was)
%       [m, n] = size(PIE.T.R.R0);                                          % MMP, 11/04/2026 (was)
%       W.dim = [0 0; m n];                                                 % MMP, 11/04/2026 (was)
%       W.I = PIE.T.I;                                                      % MMP, 11/04/2026 (was)
%       W.var1 = s1;                                                        % MMP, 11/04/2026 (was)
%       W.R.R0 = s1^k * eye(m);                                             % MMP, 11/04/2026 (was)
%       W.R.R1 = 0;                                                         % MMP, 11/04/2026 (was)
%       W.R.R2 = 0;                                                         % MMP, 11/04/2026 (was)
%                                                                           % MMP, 11/04/2026 (was)
%       PIE.T = W * PIE.T;                                                  % MMP, 11/04/2026 (was)
%       PIE.Tu = W * PIE.Tu;                                                % MMP, 11/04/2026 (was)
%       PIE.Tw = W* PIE.Tw;                                                 % MMP, 11/04/2026 (was)
        % W*X for W = s1^k I on L2: Q2, R0, R1, R2 of X times s1^k; X must  % MMP, 11/04/2026
        % have no R^n output, since W has no R^n input                      % MMP, 11/04/2026
        for op = {'T','Tu','Tw'}                                            % MMP, 11/04/2026
            X = PIE.(op{1});                                                % MMP, 11/04/2026
            if X.dim(1,1)~=0                                                % MMP, 11/04/2026
                error("The weight requires a PIE without ODE states.");     % MMP, 11/04/2026
            end                                                             % MMP, 11/04/2026
            X.Q2 = s1^k*X.Q2;                                               % MMP, 11/04/2026
            X.R.R0 = s1^k*X.R.R0;                                           % MMP, 11/04/2026
            X.R.R1 = s1^k*X.R.R1;                                           % MMP, 11/04/2026
            X.R.R2 = s1^k*X.R.R2;                                           % MMP, 11/04/2026
            PIE.(op{1}) = X;                                                % MMP, 11/04/2026
        end                                                                 % MMP, 11/04/2026
end
