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
% Add support for 2D, DJ - 07/01/2024.

% % Process the inputs.
if nargin==1
    error("Insufficient input arguments.")
elseif nargin==2
    tol = 1e-5;
end
if isa(P,'opvar2d') && isa(Z,'opvar2d')
    % Deal with 2D case.
    L = getObserver_2D(P,Z,tol);
    return
elseif ~isa(P,'opvar') || ~isa(Z,'opvar')
    error("Inputs must both be of type 'opvar', or both be of type 'opvar2d'.")
elseif isvalid(P)~=0 && isvalid(Z)~=0
    error("Input operators are not properly defined.")
end

% % Compute the actual observer gain.
L = inv(P,tol)*Z;

end