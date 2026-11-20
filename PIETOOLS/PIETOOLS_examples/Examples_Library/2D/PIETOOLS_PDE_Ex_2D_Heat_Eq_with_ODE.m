function [PDE_t, st] = PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE(GUI,params)      % MMP, 11/18/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_HEAT_EQ_WITH_ODE(GUI,PARAMS) returns the
% 2-D heat equation below, coupled to an ODE through the value x2(0,0), as
% a 'pde_struct' object, and (second output) the recommended settings ST
% of its stability test. It sets stability = 1 and stability_form = 'Q'
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
% - st:         (optional) recommended settings for the stability test      % MMP, 11/18/2026
%
% %---------------------------------------------------------------------% %
% % 2D heat equation coupled to ODE
% % ODE         x1_{t}  = (A+BK)*x1 + B*x2(s1=0,s2=0)
% % PDE         x2_{t}  = c1*x2_{(2,0)} + c2*x2_{(0,2)}
% % With BCs    x2_{(1,0)}(s1=0) = 0;     x2(s1=1) = 0;
% %             x2_{(0,1)}(s2=0) = 0;     x2(s2=1) = 0;
% %
% % Parameters A, B, K, c1 and c2 can be set.
% % Stable if A+BK is Hurwitz.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE
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
% MMP, 11/18/2026: Second output st, example-specific settings that
%                   certify the library default at fewer decision
%                   variables than the light tier (header);
% MMP, 11/19/2026: st updated to the per-example search of 11/19/2026
%                   (maintainer: no default settings): settings unchanged;
%                   the form of its tests named in the workspace.
% MMP, 11/19/2026: the settings notes of the header removed: they recorded
%                   tests and files the user does not have, not the problem;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE(GUI,params)        % MMP, 11/18/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
c1 = 1;  c2 = 1;   k = -2;
nx = 1;    ne = 1;
A = eye(nx);    B = eye(nx,ne);  K = k*eye(ne,nx);
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
ne = size(K,1);

%%% pde_var input
clear stateNameGenerator
pde_var x1 x2
x1.size = ne;
x2.vars = [s1;s2];      x2.dom = [0,1;0,1];     x2.size = ne;
PDE_t = [diff(x1,'t')==(A+B*K)*x1 +B*subs(x2,[s1;s2],[0;0]);
         diff(x2,'t')==c1*diff(x2,s1,2)+c2*diff(x2,s2,2);
         subs(diff(x2,s2),s2,0)==0;     subs(x2,s2,1)==0;
         subs(diff(x2,s1),s1,0)==0;     subs(x2,s1,1)==0];

% %%% Term-based input format
% % Initialize an ODE state component.
% PDE_t.x{1}.vars = [];
% % Initialize a 2D PDE state component.
% PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
% 
% % ODE: x1_{t} = (A+B*K)*x1
% PDE_t.x{1}.term{1}.x = 1;
% PDE_t.x{1}.term{1}.C = A+B*K;
% 
% % ODE: x1_{t} = ... + B*x2(s1=0,s2=0)
% PDE_t.x{1}.term{2}.x =2;
% PDE_t.x{1}.term{2}.loc = [0,0];
% PDE_t.x{1}.term{2}.C = B;
% 
% % PDE: x2_{t} = [c1, c2] * [x2_{s1s1}; x2_{s2s2}]
% PDE_t.x{2}.term{1}.x = 2;
% PDE_t.x{2}.term{1}.D = [2,0; 0,2];
% PDE_t.x{2}.term{1}.C = [c1*eye(ne), c2*eye(ne)];
% 
% % BC1: 0 = x2_{s2}(s1,0)
% PDE_t.BC{1}.term{1}.x = 2;
% PDE_t.BC{1}.term{1}.D = [0,1];
% PDE_t.BC{1}.term{1}.loc = [s1,0];
% % BC2: 0 = x2_{s1}(0,s2)
% PDE_t.BC{2}.term{1}.x = 2;
% PDE_t.BC{2}.term{1}.D = [1,0];
% PDE_t.BC{2}.term{1}.loc = [0,s2];
% % BC3: 0 = x2(s1,1)
% PDE_t.BC{3}.term{1}.x = 2;
% PDE_t.BC{3}.term{1}.loc = [s1,1];
% % BC4: 0 = x2(1,s2)
% PDE_t.BC{4}.term{1}.x = 2;
% PDE_t.BC{4}.term{1}.loc = [1,s2];


% Example-specific settings                                                 % MMP, 11/19/2026
% (second output only; light tier; slack lift degree at its minimum Dmin;
% Q of degree 0)
if nargout>1                                                                % MMP, 11/18/2026
    st = lpisettings_sop('light');                                          % MMP, 11/18/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
    st.sop.deg_Q{2} = 0;                                                    % MMP, 11/18/2026
end                                                                         % MMP, 11/18/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end