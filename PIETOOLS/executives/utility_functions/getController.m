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

if nargin==2
    tol = 1e-8;
end

if isvalid(P)==0 && isvalid(Z)==0
    K = Z*inv(P,tol);
else
    error("Inputs must be opvar variables");
end
end