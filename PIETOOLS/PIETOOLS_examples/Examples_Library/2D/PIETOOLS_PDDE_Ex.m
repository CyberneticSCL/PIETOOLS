function [PDE_t,st] = PIETOOLS_PDDE_Ex(GUI,params)                          % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDDE Examples
% [PDE_t,ST] = PIETOOLS_PDDE_EX(GUI,PARAMS) returns the 1-D heat equation
% with two delayed reaction terms below, the delays written as transport
% states in a second variable s2, as a 'pde_struct' object, and the
% recommended settings ST of its H2 norm test. It sets H2_norm = 1,
% H2_norm_dual = 1 and H2_norm_form = 'oco' (dual form 'cco') in the base
% workspace; a nonzero GUI prints that no GUI representation exists. The
% delays tau are set before PARAMS is evaluated, so PARAMS can change them.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% Original 1D PDE: x1_{t}(t,s1) = x1_{s1s1}(t,s1) - 20 x1(t,s1) - 4 x2(t,s1,-1) - 0.1x3(t,s1,-1) + s1(s1-\pi)w(t);
% 2D Transport Equation resulting from delay 1: x2_{t}(t,s1,s2) = 4 x2_{s2}(t,s1,s2);
%  2D Transport Equation resulting from delay 2: x3_{t}(t,s1,s2) = 3.3333 x3_{s2}(t,s1,s2);
%  (s1 in [0,pi], s2 in [-1,0]; x2(t,s1,-1) = x1(t-0.25,s1), x3(t,s1,-1) = x1(t-0.3,s1))
%  output z(t) = int_0^pi x1(t,s1) ds1;
%  BC1 0 = x1(t,0);
%  BC2 0 = x1(t,3.1416);
% BC3 0 = x2(t,s1,0) - x1(t,s1);
% BC4  0 = x3(t,s1,0) - x1(t,s1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDDE_Ex
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
% Initial coding DB - 01/20/2025
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (H2_norm_form),
%                   from the per-example search of 11/19/2026 (maintainer:
%                   no default settings);
% function [PDE_t] = PIETOOLS_PDDE_Ex(GUI,params)                           % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: tau set before the params loop, so that params can change
%   the delays; the unused pvar s t removed.

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','H2_norm = 1;');
evalin('base','H2_norm_dual = 1;');
evalin('base','H2_norm_form = ''oco'';');  % form of the H2 norm test       % MMP, 11/19/2026

% Delays: x2, x3 are transported at speeds 1/tau(1), 1/tau(2) over s2;      % MMP, 11/20/2026
% set before params, which may change them                                  % MMP, 11/20/2026
tau = [0.25;0.3];                                                           % MMP, 11/20/2026
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% Declare independent variables (time and space)                            % MMP, 11/20/2026 (was)
%pvar s t                                                                   % MMP, 11/20/2026 (was)
% Delays: x2, x3 are transported at speeds 1/tau(1), 1/tau(2) over s2       % MMP, 11/20/2026 (was)
%tau = [0.25;0.3];                                                          % MMP, 11/20/2026 (was)
% States x1(s1) and x2, x3(s1,s2); disturbance u ('input', the w above)
% and regulated output y ('output', the z above)
x1 = pde_var('state',1,s1,[0,pi]);
u = pde_var('input',1,[],[]);
y= pde_var('output',1,[],[]);
x2 = pde_var('state',1,[s1 s2],[0,pi;-1,0]);
x3 = pde_var('state',1,[s1 s2],[0,pi;-1,0]);
% Equations, BCs (Dirichlet in s1; x2 = x3 = x1 at s2 = 0) and output
PDE1=diff(x1,'t')==diff(x1,s1,2)-20*x1-4*subs(x2,s2,-1)-0.1*subs(x3,s2,-1)+(s1*(s1-pi))*u;
PDE2=diff(x2,'t')==(1/tau(1))*diff(x2,s2);
PDE3=diff(x3,'t')==(1/tau(2))*diff(x3,s2);
BC1=subs(x1,s1,0)==0;
BC2=subs(x1,s1,pi)==0;
BC3=subs(x2,s2,0)==x1;
BC4=subs(x3,s2,0)==x1;
output=y==int(x1,s1,[0,pi]);
PDE_t=[PDE1;PDE2;PDE3;BC1;BC2;BC3;BC4;output];

% Example-specific settings                                                 % MMP, 11/19/2026
% (light tier; slack lift degree at its minimum Dmin)
st = lpisettings_sop('light');                                              % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026


% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end