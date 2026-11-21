function [PDE_t,st] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_TRANSPORT_EQ_WITH_DISTURBANCE2(GUI,PARAMS)
% returns the transport equation below, with distributed disturbance
% (s-s^2)*w and integral output z, as an initialized 'pde_struct' object,
% and the recommended settings ST of its H2 norm test. It sets H2_norm =
% 1, H2_norm_dual = 1 and H2_norm_form = 'o' (dual form 'c') in the base
% workspace; a nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system using the pde_var
%               structure.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
%  PDE :        x_{t} = x_{s} + (s-s^2)w(t)
%  With BC     x(s=1) = 0
%  And output  z(t) = int(x(t,s),s,0,1)
% H2 norm computed by numerical integration of the definition: 0.1016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2
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
% Initial coding DB - 01/06/2025
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the H2 norm test; H2_norm_form = 'o' set in the
%   base workspace; the unused pvar s1 s2 removed (the PDE is in s); help
%   revised, license added.
% function [PDE_t] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);
 clear stateNameGenerator
% Initialize variables
pde_struct PDE_t;
%pvar s1 s2                                                                 % MMP, 11/21/2026 (was)

%%% Executive Function:
evalin('base','H2_norm = 1;');
evalin('base','H2_norm_dual = 1;');
evalin('base','H2_norm_form = ''o'';');  % form of the H2 norm test         % MMP, 11/21/2026

npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
% Declare the sytem equations
pde = [diff(x,t,1)==diff(x,s,1)+(s-s^2)*w;    % dynamics
                z==int(x,s,[0,1]);                     % output equation
                subs(x,s,1)==0];                            % boundary condition
PDE_t=initialize(pde);


% Example-specific settings: the veryheavy tier                             % MMP, 11/21/2026
st = lpisettings_sop('veryheavy');                                          % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end