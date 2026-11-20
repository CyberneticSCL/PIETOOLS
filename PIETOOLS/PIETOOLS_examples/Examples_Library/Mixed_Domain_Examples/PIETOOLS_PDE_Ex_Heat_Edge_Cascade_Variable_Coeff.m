function [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_HEAT_EDGE_CASCADE_VARIABLE_COEFF(GUI,PARAMS)
% returns the cascade of a 1-D and a 2-D heat equation with variable
% coefficients below, as a 'pde_struct' object, and the recommended
% settings ST of its stability test. It sets stability = 1 and
% stability_form = 'Q' (generalized storage T'Q) in the base workspace; a
% nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'r1=10;','r2=15;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings_sop structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % Cascade of a 1D and a 2D heat equation with variable coefficients:
% % the 1D state x1 is the Dirichlet datum of the 2D state x2 on the edge
% % s2 = 0.
% % PDE         x1_{t}(t,s1)    = (1+s1)*x1_{s1s1} + r1*x1;   s1 in [0,1]
% %             x2_{t}(t,s1,s2) = x2_{s1s1} + (1+s2^2)*x2_{s2s2} + r2*x2;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    x1(t,s1=0) = 0;   x1(t,s1=1) = 0;
% %             x2 = 0 at s1 = 0, 1;
% %             x2(t,s1,s2=0) = x1(t,s1);   x2(t,s1,s2=1) = 0;
% %
% % Parameters r1 and r2 can be set; default r1 = 0.9*r1*, r2 = 0.9*r2*.
% %
% % Thresholds of the blocks [finite differences, second-order central,
% % 60 interior nodes per direction]:
% %   r1* = 14.507976;   r2* = 22.205980.
% % The cascade is triangular, hence stable iff r1 < r1* and r2 < r2*
% % [inferred, not computed for the coupled system].
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff
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

% Default parameters, overridden by the assignment strings in params
r1 = 0.9*14.507976;     r2 = 0.9*22.205980;
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
x1 = pde_var(1,s1,[0,1]);
x2 = pde_var(1,[s1;s2],[0,1;0,1]);
PDE_t = [diff(x1,'t')==(1+s1)*diff(x1,s1,2)+r1*x1;
         diff(x2,'t')==diff(x2,s1,2)+(1+s2^2)*diff(x2,s2,2)+r2*x2;
         subs(x1,s1,0)==0;          subs(x1,s1,1)==0;
         subs(x2,s1,0)==0;          subs(x2,s1,1)==0;
         subs(x2,s2,0)==x1;         subs(x2,s2,1)==0];


% Example-specific settings: the light tier, with the slack integral bases
% capped at total degree theta+s <= Dmax-1
st = lpisettings_sop('light');
st.lpi.ineq.cap_inc = 'auto';

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
