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
    tol = 1e-4;
end

if isvalid(P)==0 && isvalid(Z)==0
    K = Z*inv(P,tol);
    
    %Truncate the polynomial parts on the controllers to the accuracy defined by tol.
    if ~isempty(K.Q1)
        Kp=K.Q1;
        Kpcoeff=full(Kp.coefficient);
        for i=1:size(Kpcoeff,1)
            for j=1:size(Kpcoeff,2)
                if abs(Kpcoeff(i,j))<=tol
                    Kpcoeff(i,j)=0;
                end
            end
        end
        K.Q1.coefficient=sparse(Kpcoeff);
    end
    if ~isempty(K.Q2)
        Kp=K.Q2;
        Kpcoeff=full(Kp.coefficient);
        for i=1:size(Kpcoeff,1)
            for j=1:size(Kpcoeff,2)
                if abs(Kpcoeff(i,j))<=tol
                    Kpcoeff(i,j)=0;
                end
            end
        end
        K.Q2.coefficient=sparse(Kpcoeff);
    end
    
    if ~isempty(K.R.R0)
        Kp=K.R.R0;
        Kpcoeff=full(Kp.coefficient);
        for i=1:size(Kpcoeff,1)
            for j=1:size(Kpcoeff,2)
                if abs(Kpcoeff(i,j))<=tol
                    Kpcoeff(i,j)=0;
                end
            end
        end
        K.R.R0.coefficient=sparse(Kpcoeff);
    end
    
    if ~isempty(K.R.R1)
        Kp=K.R.R1;
        Kpcoeff=full(Kp.coefficient);
        for i=1:size(Kpcoeff,1)
            for j=1:size(Kpcoeff,2)
                if abs(Kpcoeff(i,j))<=tol
                    Kpcoeff(i,j)=0;
                end
            end
        end
        K.R.R1.coefficient=sparse(Kpcoeff);
    end
    
    if ~isempty(K.R.R2)
        Kp=K.R.R2;
        Kpcoeff=full(Kp.coefficient);
        for i=1:size(Kpcoeff,1)
            for j=1:size(Kpcoeff,2)
                if abs(Kpcoeff(i,j))<=tol
                    Kpcoeff(i,j)=0;
                end
            end
        end
        K.R.R2.coefficient=sparse(Kpcoeff);
    end
else
    error("Inputs must be opvar variables");
end
end