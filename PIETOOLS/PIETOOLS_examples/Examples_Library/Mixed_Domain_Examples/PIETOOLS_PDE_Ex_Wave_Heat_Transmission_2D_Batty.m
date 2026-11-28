function [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_WAVE_HEAT_TRANSMISSION_2D_BATTY(GUI,PARAMS)
% returns the wave-heat transmission problem on two rectangles below, the
% wave reflected onto the unit square, as a 'pde_struct' object, and the
% recommended settings ST of its stability test. It sets stability = 1 and
% stability_form = 'Q' (generalized storage T'Q) in the base workspace;
% PARAMS sets nothing (no parameter), and a nonzero GUI prints that no GUI
% representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings (this example has none).
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % Wave-heat transmission on two rectangles (Batty, Paunonen, Seifert
% % 2019) [read]: wave on [-1,0]x[0,1], heat on [0,1]x[0,1]:
% % PDE         u_{tt} = u_{xx} + u_{yy};   w_{t} = w_{xx} + w_{yy};
% % With BCs    u(t,-1,y) = 0;   u = 0 at y = 0, 1;
% %             w(t,1,y) = 0;    w = 0 at y = 0, 1;
% %             u_{t}(t,0,y) = w(t,0,y);   u_{x}(t,0,y) = w_{x}(t,0,y).
% %
% % Reflect the wave onto the unit square (s1 = -x, s2 = y) and let
% % v = u_{t}; v is 1st order in s1 so that its trace at s1 = 0 exists:
% % PDE         u_{t} = v;   v_{t} = u_{s1s1} + u_{s2s2};
% %             w_{t} = w_{s1s1} + w_{s2s2};       (s1,s2) in [0,1]^2
% % With BCs    u(t,1,s2) = 0;   v(t,1,s2) = 0;   w(t,1,s2) = 0;
% %             u_{s1}(t,0,s2) = -w_{s1}(t,0,s2);   w(t,0,s2) = v(t,0,s2);
% %             u = v = w = 0 at s2 = 0, 1.
% %
% % Known result [read]: the energy is non-increasing and decays as
% % o(t^{-2/3}) for classical solutions, a sharp rate; the semigroup is
% % not exponentially stable (rays parallel to the interface). No
% % parameter threshold.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty
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
u = pde_var([s1;s2],[0,1;0,1],[2;2]);
v = pde_var([s1;s2],[0,1;0,1],[1;2]);           % 1st order in s1, 2nd in s2
w = pde_var([s1;s2],[0,1;0,1],[2;2]);
PDE_t = [diff(u,'t')==v;
         diff(v,'t')==diff(u,s1,2)+diff(u,s2,2);
         diff(w,'t')==diff(w,s1,2)+diff(w,s2,2);
         subs(u,s1,1)==0;           subs(diff(u,s1),s1,0)==-subs(diff(w,s1),s1,0);
         subs(v,s1,1)==0;
         subs(w,s1,1)==0;           subs(w,s1,0)==subs(v,s1,0);
         subs([u;v;w],s2,0)==0;     subs([u;v;w],s2,1)==0];


% Example-specific settings: the light tier
st = lpisettings('light');

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{batty2019optimal,
%   title={Optimal energy decay for the wave-heat system on a rectangular domain},
%   author={Batty, Charles J K and Paunonen, Lassi and Seifert, David},
%   journal={SIAM Journal on Mathematical Analysis},
%   volume={51},
%   number={2},
%   pages={808--819},
%   year={2019},
%   doi={10.1137/18M1195796},
%   note={arXiv:1803.05405}
% }
