function [PDE_t,st] = PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_3D_ADVECTION_DIFFUSION_EQ(GUI,PARAMS)
% returns the 3-D advection-diffusion-reaction equation below, with
% Dirichlet BCs, as a 'pde_struct' object, and the recommended settings ST
% of its stability test. It sets stability = 1 and stability_form = 'P'
% (storage <T x, P T x>, P >= eppos I) in the base workspace; a nonzero
% GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','C=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - st:         recommended settings (lpisettings_sop structure) for the
%               stability test of this example.
%
% %---------------------------------------------------------------------% %
% % 3D advection-diffusion-reaction equation (3D version of the equation of
% % Holmes, 1994, with a reaction term):
% % PDE         x_{t}   = C*(x_{s1s1} + x_{s2s2} + x_{s3s3})
% %                       + b1*x_{s1} + b2*x_{s2} + b3*x_{s3} + lam*x
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;   x(s3=0) = 0;
% %             x(s1=1) = 0;   x(s2=1) = 0;   x(s3=1) = 0;
% %
% % Parameters C, b1, b2, b3, lam can be set.
% %
% % Substituting x = exp(-(b1*s1+b2*s2+b3*s3)/(2*C))*y gives the reaction-
% % diffusion equation y_{t} = C*(y_{s1s1}+y_{s2s2}+y_{s3s3})
% % + (lam - (b1^2+b2^2+b3^2)/(4*C))*y, with y=0 on the boundary, with
% % eigenvalues lam - (b1^2+b2^2+b3^2)/(4*C) - C*pi^2*(k1^2+k2^2+k3^2).
% % The system is therefore exponentially stable if and only if
% %   lam < 3*pi^2*C + (b1^2+b2^2+b3^2)/(4*C),
% % i.e. lam < 30.4838 for the default C = 1, b = (0.5, 1, 1.5). The
% % advection stabilizes the system: without it, lam < 3*pi^2 = 29.6088.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 11/02/2026

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2 s3

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
C = 1;  b1 = 0.5;  b2 = 1;  b3 = 1.5;  lam = 1;
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
PDE_t = [diff(x,'t')==C*(diff(x,s1,2)+diff(x,s2,2)+diff(x,s3,2)) ...
                      + b1*diff(x,s1)+b2*diff(x,s2)+b3*diff(x,s3) + lam*x;
         subs(x,s1,0)==0;   subs(x,s1,1)==0;
         subs(x,s2,0)==0;   subs(x,s2,1)==0;
         subs(x,s3,0)==0;   subs(x,s3,1)==0];


% Example-specific settings: the stripped tier
st = lpisettings_sop('stripped');

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, "Partial differential
% equations in ecology: Spatial interactions and population dynamics,"
% Ecology, vol. 75, pp. 17-29, 1994.