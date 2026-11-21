function [PDE_t,st] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_UNSTABLE_REACTIONDIFFUSION(GUI,PARAMS)
% returns the unstable reaction-diffusion equation below, with
% disturbance w, regulated output z and measured output y, as a
% 'pde_struct' object, and the recommended settings ST of its H2 optimal
% estimator synthesis. It sets H2_estimator = 1 in the base workspace;
% PARAMS is not used. With GUI = 1 it loads a GUI save file into a new
% PIETOOLS_PDE_GUI window.
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
%  PDE :        x_{t} = 3x+(s^2+0.2)x_{ss} -(s^2)/2w(t)
%  With BC     x(s=0) = 0, x_{s}(s=1)=0
%  And outputs  z(t) = int(x(t,s),s,0,1)
%                       y(t)=x(s=1)+w
% % 
% % No parameter can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion
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
% Initial coding DJ - 08/18/2022
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function       [PDE_t] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion(GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).
% MMP, 11/20/2026: st.sop.deg_Q removed: only the stability executive reads
%   it, and this example names the H2 estimator synthesis.

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','H2_estimator = 1;');

% Construct the PDE.
pvar s t;
x = pde_var('state',1,s,[0,1]); 
z = pde_var('out',1); 
y=pde_var('sense',1); 
w=pde_var('in',1); 
PDE_t = [diff(x,t)==3*x+(s^2+0.2)*diff(x,s,2)-0.5*s^2*w; subs(x,s,0)==0; subs(diff(x,s),s,1)==0;
         z==int(x,s,[0,1]); y==subs(x,s,1)+w];
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/19/2026
% (heavy tier; 1-D storage with upper-integral s-length 1 and theta-lengths
% 3, 2, 1 (multiplier, upper, lower integral); free operators of degree 2;
% no slack weight raise of the synthesis LPI; Q of degree 6 in the          % MMP, 11/20/2026 (was)
% multiplier and 2 in the integral kernels)                                 % MMP, 11/20/2026 (was)
% no slack weight raise of the synthesis LPI)                               % MMP, 11/20/2026
st = lpisettings_sop('heavy');                                              % MMP, 11/19/2026
st.lpi.pos.len_com{1} = 1;                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [3;2;1];                                            % MMP, 11/19/2026
st.lpi.var.deg{1} = [2 2 2];                                                % MMP, 11/19/2026
st.sop.dw_P = 0;                                                            % MMP, 11/19/2026
%st.sop.deg_Q{1} = [6 2 2];                                                 % MMP, 11/20/2026 (was)

end
% @article{peet2019discussion,
%   title={Discussion paper: A new mathematical framework for representation and analysis of coupled pdes},
%   author={Peet, Matthew M and Shivakumar, Sachin and Das, Amritam and Weiland, Seip},
%   journal={IFAC-PapersOnLine},
%   volume={52},
%   number={2},
%   pages={132--137},
%   year={2019},
%   publisher={Elsevier}
% }