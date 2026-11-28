function [PDE_t, st] = PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_HEAT_LAPLACE_BELTRAMI_EDGE(GUI,PARAMS)
% returns the 2-D heat equation with a dynamic (Laplace-Beltrami) boundary
% condition on the edge s2 = 1 below, the edge state v in s1, as a
% 'pde_struct' object, and (second output) the recommended settings ST of
% its stability test. It sets stability = 1 and stability_form = 'Q'
% (generalized storage T'Q) in the base workspace; a nonzero GUI prints
% that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=0.5;','sig=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) recommended settings (lpisettings
%               structure) for the stability test of this example.
%
% %---------------------------------------------------------------------% %
% % 2D heat equation with a dynamic boundary condition of reactive-
% % diffusive (Laplace-Beltrami) type on the edge s2 = 1 (box version of
% % Goldstein 2006; Vazquez, Vitillaro 2011) [derived]:
% % PDE         u_{t}(t,s1,s2) = u_{s1s1} + u_{s2s2} + lam*u;
% %                                               (s1,s2) in [0,1]^2
% %             v_{t}(t,s1)    = sig*v_{s1s1} - u_{s2}(t,s1,s2=1) + lam*v;
% %                                               s1 in [0,1]
% % With BCs    u_{s1} = 0 at s1 = 0, 1;   v_{s1} = 0 at s1 = 0, 1;
% %             u(t,s1,s2=0) = 0;          u(t,s1,s2=1) = v(t,s1).
% % The edge state v is the Dirichlet datum of u on s2 = 1, and the flux
% % of u drives v.
% %
% % Parameters lam and sig can be set; default sig = 1, lam = 0.9*lam*.
% %
% % Threshold [derived from the modes, exact]: in the cosine modes in s1,
% % mode n = 0 gives u = sin(k*s2), v = sin(k), with k*tan(k) = 1; modes
% % n >= 1 subtract (n*pi)^2 from every eigenvalue (sig = 1). Hence
% %   stable iff lam < lam* = k0^2 = 0.740174,   k0 = 0.860334.
% % No threshold is published for a box.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge
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
% Initial coding MMP - 11/17/2026

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test

% Default parameters (lam = 0.9*lam*), overridden by the strings in params
sig = 1;    lam = 0.9*0.740174;
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
v = pde_var(s1,[0,1]);
u = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(v,'t')==sig*diff(v,s1,2)-subs(diff(u,s2),s2,1)+lam*v;
         diff(u,'t')==diff(u,s1,2)+diff(u,s2,2)+lam*u;
         subs(diff(v,s1),s1,0)==0;  subs(diff(v,s1),s1,1)==0;
         subs(diff(u,s1),s1,0)==0;  subs(diff(u,s1),s1,1)==0;
         subs(u,s2,0)==0;           subs(u,s2,1)==v];


% Example-specific settings (second output only): the light tier, with
if nargout>1
    st = lpisettings('light');
    st.lpi.ineq.dD = 0;             % slack lift degree at its minimum Dmin
    st.lpi.ineq.cap_inc = 'auto';   % slack integral bases: theta+s <= Dmax-1
end

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{goldstein2006derivation,
%   title={Derivation and physical interpretation of general boundary conditions},
%   author={Goldstein, Gis{\`e}le Ruiz},
%   journal={Advances in Differential Equations},
%   volume={11},
%   number={4},
%   pages={457--480},
%   year={2006}
% }
% @article{vazquez2011heat,
%   title={Heat equation with dynamical boundary conditions of reactive-diffusive type},
%   author={V{\'a}zquez, Juan Luis and Vitillaro, Enzo},
%   journal={Journal of Differential Equations},
%   volume={250},
%   pages={2143--2161},
%   year={2011},
%   doi={10.1016/j.jde.2010.12.012}
% }
% (Goldstein 2006: one repository lists pages 419-456; the DOI of
% Vazquez, Vitillaro 2011 is unverified.)
