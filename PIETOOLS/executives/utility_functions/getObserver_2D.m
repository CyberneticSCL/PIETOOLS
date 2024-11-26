function Lop = getObserver_2D(Pop,Zop,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getObserver.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function returns oberver gains L = P^{-1}Z as an opvar2d object given
% Lyapunov opvar variable P and opvar2d variable Z used in the LPI for Hinf
% optimal observer LPI

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following inputs must be passed:
%
% Pop, Zop - opvar2d objects with matching inner dimensions
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 07/01/2024

% % % Process the inputs;
if nargin==1
    error("Insufficient input arguments.")
elseif nargin==2
    tol = 1e-5;
end
if ~isa(Pop,'opvar2d') || ~isa(Zop,'opvar2d')
    error("Inputs must both be of type 'opvar2d'.")
end

% % If the operator is separable, take the inverse and multiply.
if is_separable(Pop)
    Lop = clean_opvar(inv(Pop,tol),tol)*Zop;
else
    % Otherwise, try computing P\Z 
    % (works quite well, but is very expensive...)
    Lop = mldivide(Pop,Zop,2*ones(4,2),tol);
end
end