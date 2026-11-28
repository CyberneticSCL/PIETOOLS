function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_DELAYED_HEAT_EQ_FRIDMAN(GUI,PARAMS)
% returns the 2-D heat equation with an interior delay below, the delay
% written as a transport state x2(s1,s2,s3), as a 'pde_struct' object, and
% the recommended settings ST of its stability test. It sets stability = 1
% and stability_form = 'Q' (generalized storage T'Q) in the base
% workspace; a nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'r=2.9;','tau=0.5;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % 2D heat equation with an interior delay (2D extension of Fridman,
% % Orlov 2009, Example 1):
% % PDE         x_{t}(t,s1,s2) = x_{s1s1} + x_{s2s2} + r*x(t,s1,s2)
% %                                - x(t-tau,s1,s2);    (s1,s2) in [0,pi]^2
% % With BCs    x = 0 at s1 = 0, s1 = pi, s2 = 0, s2 = pi.
% %
% % Let x1 = x, and x2(t,s1,s2,s3) = x1(t-s3,s1,s2), s3 in [0,tau]:
% % PDE         x1_{t} = x1_{s1s1} + x1_{s2s2} + r*x1 - x2(t,s1,s2,s3=tau);
% %             x2_{t} = -x2_{s3};
% % With BCs    x1 = x2 = 0 at s1 = 0, pi and at s2 = 0, pi;
% %             x2(t,s1,s2,s3=0) = x1(t,s1,s2).
% % The 2D state x1 is the datum of the 3D state x2 on the face s3 = 0.
% %
% % Parameters r and tau can be set; default r = 2.9, tau = 0.9*tau*.
% %
% % Threshold [derived from the modes, exact]: the sine mode (m,n) obeys
% % y' = a*y - y(t-tau), a = r-m^2-n^2, with delay margin
% % arccos(a)/sqrt(1-a^2) for a in (-1,1), and stable for every delay for
% % a <= -1. For r = 2.9 the minimum is attained at (m,n) = (1,1), a = 0.9:
% %   stable iff tau < tau* = 1.034726,
% % the threshold of the 1D version with r = 1.9, which
% % Fridman, Orlov 2009, Remark 4, show to be tight [read].
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman
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
pvar s1 s2 s3

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test

% Default parameters (tau = 0.9*tau*), overridden by the strings in params
r = 2.9;    tau = 0.9*1.034726;
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
x1 = pde_var([s1;s2],[0,pi;0,pi]);
x2 = pde_var([s1;s2;s3],[0,pi;0,pi;0,tau],[2;2;1]);   % 2nd order in s1, s2; 1st in s3
PDE_t = [diff(x1,'t')==diff(x1,s1,2)+diff(x1,s2,2)+r*x1-subs(x2,s3,tau);
         diff(x2,'t')==-diff(x2,s3);
         subs([x1;x2],s1,0)==0;     subs([x1;x2],s1,pi)==0;
         subs([x1;x2],s2,0)==0;     subs([x1;x2],s2,pi)==0;
         subs(x2,s3,0)==x1];


% Example-specific settings: the light tier
st = lpisettings('light');

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{fridman2009exponential,
%   title={Exponential stability of linear distributed parameter systems with time-varying delays},
%   author={Fridman, Emilia and Orlov, Yury},
%   journal={Automatica},
%   volume={45},
%   number={1},
%   pages={194--201},
%   year={2009},
%   doi={10.1016/j.automatica.2008.06.006}
% }
