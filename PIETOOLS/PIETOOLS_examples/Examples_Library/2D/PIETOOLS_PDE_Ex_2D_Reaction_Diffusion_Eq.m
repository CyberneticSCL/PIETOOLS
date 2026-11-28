function [PDE_t, st] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_REACTION_DIFFUSION_EQ(GUI,PARAMS)
% returns the 2-D advection-diffusion equation below, with Dirichlet BCs,
% as a 'pde_struct' object, and (second output) the recommended settings
% ST of its stability test. It sets stability = 1 and stability_form = 'P'
% (storage <T x, P T x>, P >= eppos I) in the base workspace; a nonzero
% GUI prints that no GUI representation exists.
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
% - st:         (optional) recommended settings for the stability test      % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % 2D advection-diffusion equation (Holmes, 1994):
% % PDE         x_{t}   = C*(x_{s1s1} + x_{s2s2}) + b1*x_{s1} + b2*x_{s2}
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% %             x(s1=1) = 0;   x(s2=1) = 0;
% %
% % Parameters C, b1, and b2 can be set.
% % Exponentially stable for C > 0 (eigenvalues
% %   -C*pi^2*(m^2+n^2) - (b1^2+b2^2)/(4*C), m, n >= 1).
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq
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
% MMP, 11/21/2026: Changes since PIETOOLS 2025: optional second output st,
%   the recommended settings of the stability test; stability_form = 'P'
%   set in the base workspace; help revised (OUTPUT lists no PDE_b, which
%   is not returned; stability statement added), license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params)   % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/21/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
C = 1;  b1 = 0.5;  b2 = 2;  ne = 1;
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
x = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(x,'t')==C*(diff(x,s1,2)+diff(x,s2,2)) +b1*diff(x,s1)+b2*diff(x,s2);
         subs(x,s1,0)==0;   subs(x,s1,1)==0;
         subs(x,s2,0)==0;   subs(x,s2,1)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
% 
% % PDE: x_{t} = [C, C] * [x_{s1s1}; x_{s2s2}];
% PDE_t.x{1}.term{1}.D = [2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [C*eye(ne), C*eye(ne)];
% 
% % PDE: x_{t} = ... + [b1, b2] * [x_{s1}; x_{s2}];
% PDE_t.x{1}.term{2}.D = [1,0; 0,1];
% PDE_t.x{1}.term{2}.C = [b1*eye(ne), b2*eye(ne)];
% 
% % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
% PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
% PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];


% Example-specific settings                                                 % MMP, 11/21/2026
% (second output only; stripped tier; slack lift degree at its minimum
% Dmin; storage on the 2-D space with theta-lengths [0 1] for the
% multiplier (none in s1) and [1 1] for the integral; 2-D slack with the
% plain term only, psatz 'none')
if nargout>1                                                                % MMP, 11/21/2026
    st = lpisettings('stripped');                                           % MMP, 11/21/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/21/2026
    st.lpi.pos.len_sos{2} = [0 1;1 1;0 0];                                  % MMP, 11/21/2026
    st.lpi.ineq.psatz{2} = 'none';                                          % MMP, 11/21/2026
end                                                                         % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.