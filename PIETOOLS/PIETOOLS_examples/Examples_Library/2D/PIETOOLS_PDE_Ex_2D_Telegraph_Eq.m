function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Telegraph_Eq(GUI,params)           % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_TELEGRAPH_EQ(GUI,PARAMS) returns the 2-D
% telegraph equation below, in the states x1 = x, x2 = x_t, with Dirichlet
% BCs, as a 'pde_struct' object, and the recommended settings ST of its
% stability test. It sets stability = 1 and stability_form = 'Q'
% (generalized storage T'Q) in the base workspace; a nonzero GUI prints
% that no GUI representation exists.
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
% % % Telegraph equation (Holmes, 1994):
% % % PDE         x_{t}   = -1/(2*lam) * x_{tt} + (C^2/(2*lam))*(x_{s1s1} + x_{s2s2})
% % % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% % %             x(s1=1) = 0;   x(s2=1) = 0;
% % %
% % % Use states x1 = x, x2 = x_{t}, then:
% % %
% % % PDE         x1_{t}   = x2
% % %             x2_{t}   = -2*lam * x2 + C^2*(x1_{s1s1} + x1_{s2s2})
% % % With BCs    x1(s1=0) = 0;   x1(s2=0) = 0;
% % %             x1(s1=1) = 0;   x1(s2=1) = 0;
% % %             x2(s1=0) = 0;   x2(s2=0) = 0;
% % %             x2(s1=1) = 0;   x2(s2=1) = 0;
% %
% % Parameters C, lam and the state size ne can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_Telegraph_Eq
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
%   returned; coefficient C^2/(2*lam) corrected), license added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_Telegraph_Eq(GUI,params)            % MMP, 11/21/2026 (was)

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
C = 1;  lam = 1;    ne = 1;
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
x1 = pde_var(ne,[s1;s2],[0,1;0,1]);
x2 = pde_var(ne,[s1;s2],[0,1;0,1],[2;2]); % x2 must be 2nd order differentiable wrt both s1 and s2
PDE_t = [diff(x1,'t')==x2;
         diff(x2,'t')==-2*lam*x2 +C^2*(diff(x1,s1,2)+diff(x1,s2,2));
         subs([x1;x2],s1,0)==0;     subs([x1;x2],s1,1)==0;
         subs([x1;x2],s2,0)==0;     subs([x1;x2],s2,1)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
% PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
% PDE_t.x{2}.diff = [2,2];    % x2 must be 2nd order differentiable wrt both s1 and s2
% 
% % PDE: x1_{t} = x2;
% PDE_t.x{1}.term{1}.x = 2;
% 
% % PDE: x2_{t} = -2*lam*x2
% PDE_t.x{2}.term{1}.x = 2;
% PDE_t.x{2}.term{1}.C = -2*lam;
% 
% % PDE: x2_{t} = ... + [C^2, C^2] * [x1_{s1s1}; x1_{s2s2}];
% PDE_t.x{2}.term{2}.x = 1;
% PDE_t.x{2}.term{2}.D = [2,0; 0,2];
% PDE_t.x{2}.term{2}.C = [C^2 * eye(ne), C^2 * eye(ne)];
% 
% % BC1: 0 = x1(s1,0)                    % BC3: 0 = x1(s1,1)
% PDE_t.BC{1}.term{1}.x = 1;             PDE_t.BC{3}.term{1}.x = 1;
% PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% % BC2: 0 = x1(0,s2)                    % BC4: 0 = x1(1,s2)
% PDE_t.BC{2}.term{1}.x = 1;             PDE_t.BC{4}.term{1}.x = 1;
% PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];
% 
% % BC5: 0 = x1(s1,0)                    % BC7: 0 = x1(s1,1)
% PDE_t.BC{5}.term{1}.x = 2;             PDE_t.BC{7}.term{1}.x = 2;
% PDE_t.BC{5}.term{1}.loc = [s1,0];      PDE_t.BC{7}.term{1}.loc = [s1,1];
% % BC6: 0 = x1(0,s2)                    % BC8: 0 = x1(1,s2)
% PDE_t.BC{6}.term{1}.x = 2;             PDE_t.BC{8}.term{1}.x = 2;
% PDE_t.BC{6}.term{1}.loc = [0,s2];      PDE_t.BC{8}.term{1}.loc = [1,s2];


% Example-specific settings: the light tier                                 % MMP, 11/21/2026
st = lpisettings('light');                                                  % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.