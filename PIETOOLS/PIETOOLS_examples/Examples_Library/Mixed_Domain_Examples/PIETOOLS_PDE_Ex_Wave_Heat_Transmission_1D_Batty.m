function [PDE_t,st] = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_WAVE_HEAT_TRANSMISSION_1D_BATTY(GUI,PARAMS)
% returns the 1-D wave-heat transmission problem below, the wave
% reflected onto [0,1] in Riemann form, as a 'pde_struct' object, and the
% recommended settings ST of its stability test. It sets stability = 1 and
% stability_form = 'P' (storage <T x, P T x>, P >= eppos I) in the base
% workspace; PARAMS sets nothing (no parameter), and a nonzero GUI prints
% that no GUI representation exists.
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
% - st:         recommended settings (lpisettings_sop structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % Wave-heat transmission in 1D (Batty, Paunonen, Seifert 2016, Dirichlet
% % variant of Sec. 5; Zhang, Zuazua 2004) [read]:
% % PDE         u_{tt} = u_{xx}  on x in [-1,0];   w_{t} = w_{xx}  on x in [0,1]
% % With BCs    u(t,-1) = 0;   w(t,1) = 0;
% %             u_{t}(t,0) = w(t,0);   u_{x}(t,0) = w_{x}(t,0).
% %
% % Reflect the wave onto [0,1] (s1 = -x) and let e = u_{s1}, v = u_{t}:
% % PDE         e_{t} = v_{s1};   v_{t} = e_{s1};   w_{t} = w_{s1s1};
% %                                               s1 in [0,1]
% % With BCs    v(t,s1=1) = 0;   e(t,s1=0) = -w_{s1}(t,0);
% %             w(t,s1=1) = 0;   w(t,s1=0) = v(t,s1=0).
% %
% % Known result [read]: the energy is non-increasing and decays as
% % o(t^{-4}) for classical solutions, a sharp rate; the semigroup is not
% % exponentially stable (spectrum in the open left half-plane,
% % approaching the imaginary axis). No parameter threshold.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty
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
pvar s1

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test

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
e = pde_var(s1,[0,1],1);
v = pde_var(s1,[0,1],1);
w = pde_var(s1,[0,1]);
PDE_t = [diff(e,'t')==diff(v,s1);
         diff(v,'t')==diff(e,s1);
         diff(w,'t')==diff(w,s1,2);
         subs(v,s1,1)==0;           subs(e,s1,0)==-subs(diff(w,s1),s1,0);
         subs(w,s1,1)==0;           subs(w,s1,0)==subs(v,s1,0)];


% Example-specific settings: the stripped tier, with the 1-D storage from
% the lower-integral basis operator only (multiplier theta-length 0,
% upper-integral s-length 0), theta-length 1; slack lift degree at its
% minimum Dmin
st = lpisettings_sop('stripped');
st.lpi.pos.len_com{1} = 0;
st.lpi.pos.len_sos{1} = [0;2;1];
st.lpi.ineq.dD = 0;

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{batty2016optimal,
%   title={Optimal energy decay in a one-dimensional coupled wave-heat system},
%   author={Batty, Charles J K and Paunonen, Lassi and Seifert, David},
%   journal={Journal of Evolution Equations},
%   volume={16},
%   pages={649--664},
%   year={2016},
%   note={arXiv:1505.04951; DOI 10.1007/s00028-015-0316-y unverified}
% }
% @article{zhang2004polynomial,
%   title={Polynomial decay and control of a 1-d hyperbolic-parabolic coupled system},
%   author={Zhang, Xu and Zuazua, Enrique},
%   journal={Journal of Differential Equations},
%   volume={204},
%   pages={380--438},
%   year={2004},
%   doi={10.1016/j.jde.2004.02.004}
% }
