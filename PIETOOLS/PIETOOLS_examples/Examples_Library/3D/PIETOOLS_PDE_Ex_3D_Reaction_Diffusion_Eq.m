function [PDE_t] = PIETOOLS_PDE_Ex_3D_Reaction_Diffusion_Eq(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - PDE_b:      PDE structure defining the example system in the batch
%               format. NOT SUPPORTED FOR ND PDEs
%
% %---------------------------------------------------------------------% %
% % 3D reaction-diffusion equation:
% % PDE         x_{t}   = C*(x_{s1s1} + x_{s2s2} + x_{s3s3}) + lam*x
% % With BCs    x(s1=0) = 0;        x(s1=1) = 0;
% %             x(s2=0) = 0;        x_{s2}(s2=1) = 0;
% %             x_{s3}(s3=0) = 0;   x(s3=1) = 0;
% %
% % Parameters C, lam can be set.
% %
% % The eigenfunctions are products of sin(k1*pi*s1), sin((k2-1/2)*pi*s2),
% % cos((k3-1/2)*pi*s3), with eigenvalues
% %   lam - C*pi^2*(k1^2 + (k2-1/2)^2 + (k3-1/2)^2),   k1,k2,k3 = 1,2,...
% % so the largest is lam - C*(pi^2 + pi^2/4 + pi^2/4) = lam - 3*pi^2*C/2.
% % The system is exponentially stable if and only if
% %   lam < 3*pi^2*C/2,      i.e. lam < 14.8044 for C = 1.
% %---------------------------------------------------------------------% %
%
% DJ, 11/02/2026: Initial coding;

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2 s3

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

C = 1;  lam = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input
clear stateNameGenerator
x = pde_var([s1;s2;s3],[0,1;0,1;0,1]);
PDE_t = [diff(x,'t')==C*(diff(x,s1,2)+diff(x,s2,2)+diff(x,s3,2)) + lam*x;
         subs(x,s1,0)==0;               subs(x,s1,1)==0;
         subs(x,s2,0)==0;               subs(diff(x,s2),s2,1)==0;
         subs(diff(x,s3),s3,0)==0;      subs(x,s3,1)==0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end