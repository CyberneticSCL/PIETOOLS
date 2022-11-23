function [L] = getObserver(P,Z,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getObserver.m     PIETOOLS 2021a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns oberver gains L = P^{-1}Z as an opvar object given
% Lyapunov opvar variable P and opvar variable Z used in the LPI for hinf
% optimal observer LPI

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
L = inv(P,tol)*Z;
else
    error("Inputs must be opvar variables");
end
end