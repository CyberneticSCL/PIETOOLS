function [PDE_t,st] = PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_3D_COUPLED_HEAT_EQ(GUI,PARAMS) returns the
% two coupled 3-D heat equations below, with Dirichlet BCs, as a
% 'pde_struct' object, and the recommended settings ST of its stability
% test. It sets stability = 1 and stability_form = 'P' (storage
% <T x, P T x>, P >= eppos I) in the base workspace; a nonzero GUI prints
% that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','C1=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - st:         recommended settings (lpisettings structure) for the
%               stability test of this example.
%
% %---------------------------------------------------------------------% %
% % 3D coupled heat equations, with symmetric coupling lam and different
% % diffusion coefficients C1, C2:
% % PDE         x1_{t}  = C1*(x1_{s1s1} + x1_{s2s2} + x1_{s3s3}) + lam*x2
% %             x2_{t}  = C2*(x2_{s1s1} + x2_{s2s2} + x2_{s3s3}) + lam*x1
% % With BCs    xi(s1=0) = 0;   xi(s2=0) = 0;   xi(s3=0) = 0;
% %             xi(s1=1) = 0;   xi(s2=1) = 0;   xi(s3=1) = 0;    i = 1,2.
% %
% % Parameters C1, C2, lam can be set.
% %
% % For each eigenfunction of the Dirichlet Laplacian, with eigenvalue
% % -mu, mu = pi^2*(k1^2+k2^2+k3^2) >= 3*pi^2, the coefficients of x1 and x2
% % evolve by the matrix [-C1*mu, lam; lam, -C2*mu], which is Hurwitz if and
% % only if C1*C2*mu^2 > lam^2. The system is therefore exponentially stable
% % if and only if
% %   |lam| < 3*pi^2*sqrt(C1*C2),
% % i.e. |lam| < 16.2174 for the defaults C1 = 1, C2 = 0.3.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq
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
C1 = 1;  C2 = 0.3;  lam = 1;
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
x1 = pde_var([s1;s2;s3],[0,1;0,1;0,1]);
x2 = pde_var([s1;s2;s3],[0,1;0,1;0,1]);
PDE_t = [diff(x1,'t')==C1*(diff(x1,s1,2)+diff(x1,s2,2)+diff(x1,s3,2)) + lam*x2;
         diff(x2,'t')==C2*(diff(x2,s1,2)+diff(x2,s2,2)+diff(x2,s3,2)) + lam*x1;
         subs(x1,s1,0)==0;   subs(x1,s1,1)==0;
         subs(x1,s2,0)==0;   subs(x1,s2,1)==0;
         subs(x1,s3,0)==0;   subs(x1,s3,1)==0;
         subs(x2,s1,0)==0;   subs(x2,s1,1)==0;
         subs(x2,s2,0)==0;   subs(x2,s2,1)==0;
         subs(x2,s3,0)==0;   subs(x2,s3,1)==0];


% Example-specific settings: the stripped tier
st = lpisettings('stripped');

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end