function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params)             % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_KISS_MODEL(GUI,PARAMS) returns the 2-D
% KISS model below, a reaction-diffusion equation with Dirichlet BCs, as a
% 'pde_struct' object, and the recommended settings ST of its stability
% test. It sets stability = 1 and stability_form = 'Q' (generalized
% storage T'Q) in the base workspace; a nonzero GUI prints that no GUI
% representation exists.
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
% - st:         example-specific settings (lpisettings structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % KISS Model from Holmes (1994):
% % PDE         x_{t}   = lam*x + c1*x_{s1s1} + c2*x_{s2s2}
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% %             x(s1=1) = 0;   x(s2=1) = 0;
% %
% % Parameters c1, c2, lam and the state size ne can be set.
% % Stable for lam <= 2*pi^2 when c1=c2=1.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_KISS_Model
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
% Initial coding DJ - 08/20/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the stability test; stability_form = 'Q' set in
%   the base workspace; help revised (OUTPUT lists no PDE_b, which is not
%   returned), license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params)              % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/21/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
c1 = 1; c2 = 1; lam = 19;   ne = 1;
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
x = pde_var(ne,[s1;s2],[0,1;0,1]);
PDE_t = [diff(x,'t')==lam*x+c1*diff(x,s1,2)+c2*diff(x,s2,2);
         subs(x,s1,0)==0;   subs(x,s1,1)==0;
         subs(x,s2,0)==0;   subs(x,s2,1)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
% 
% % PDE: x_{t} = [lam, c1, c2] * [x; x_{s1s1}; x_{s2s2}]
% PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [lam*eye(ne), c1*eye(ne), c2*eye(ne)];
% 
% % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
% PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
% PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];


% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; Q of degree 0 in the multiplier, 1 in the integrals)
st = lpisettings('stripped');                                               % MMP, 11/21/2026
st.sop.deg_Q{2} = [0 1 1];                                                  % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.