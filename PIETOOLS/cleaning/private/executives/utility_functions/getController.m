function [K] = getController(P,Z,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getController.m     PIETOOLS 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns controller gains K = ZP^{-1} as an opvar object given
% Lyapunov opvar variable P and opvar variable Z used in the LPI for hinf
% optimal controller LPI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be defined passed:
%
% P, Z - opvar variable with matching inner dimensions
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 5/20/2021
% MMP, 10/09/2026: tol is no longer passed to inv. There it cut the Neumann
%                  series and zeroed the coefficients of the inverse at 1e-4,
%                  which left |K P - Z|/|Z| at 8e-5 on the test plant of
%                  sopvar/Testfolder/Test_copvar_inv; the inverse now runs at
%                  its own default (1e-8) and tol only truncates K, as the
%                  header says. The container route is getController_sop.

if nargin==2
    tol = 1e-4;
end

if isvalid(P)==0 && isvalid(Z)==0
%   K = Z*inv(P,tol);                                                       % MMP, 10/09/2026 (was)
    K = Z*inv(P);                                                           % MMP, 10/09/2026
    
    %Truncate the polynomial parts on the controllers to the accuracy defined by tol.
    K = clean_opvar(K,tol);
else
    error("Inputs must be opvar variables");
end
end