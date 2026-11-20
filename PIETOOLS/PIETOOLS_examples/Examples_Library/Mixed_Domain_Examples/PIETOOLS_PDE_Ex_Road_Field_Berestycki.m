function [PDE_t,st] = PIETOOLS_PDE_Ex_Road_Field_Berestycki(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_ROAD_FIELD_BERESTYCKI(GUI,PARAMS) returns the
% linearized road-field model below, a road state r(s1) on the edge s2 = 0
% of a field state f(s1,s2), as a 'pde_struct' object, and the recommended
% settings ST of its stability test. It sets stability = 1 and
% stability_form = 'P' (storage <T x, P T x>, P >= eppos I) in the base
% workspace; a nonzero GUI prints that no GUI representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=10;','Dr=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings_sop structure) for the
%               test this file names.
%
% %---------------------------------------------------------------------% %
% % Road-field model (Berestycki, Roquejoffre, Rossi 2013, eq. (1.1),
% % linearized, on a box [derived]): a road r(s1) on the edge s2 = 0 of a
% % field f(s1,s2):
% % PDE         r_{t}(t,s1)    = Dr*r_{s1s1} + nu*f(t,s1,s2=0) - mu*r;
% %                                               s1 in [0,1]
% %             f_{t}(t,s1,s2) = d*(f_{s1s1} + f_{s2s2}) + lam*f;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    r = 0 at s1 = 0, 1;   f = 0 at s1 = 0, 1 and at s2 = 1;
% %             -d*f_{s2}(t,s1,0) = mu*r(t,s1) - nu*f(t,s1,0);
% % lam is the linearized growth rate f'(0) of the field.
% %
% % Parameters lam, Dr, d, mu and nu can be set; default Dr = 2, d = 1,
% % mu = nu = 1, lam = 0.9*lam*.
% %
% % The paper gives a spreading-speed threshold (Dr = 2*d) [read], not a
% % stability threshold. Threshold [derived, finite differences]: sine
% % modes in s1, finite differences in s2 with a ghost node for the Robin
% % condition, modes n = 1..3, default Dr, d, mu, nu:
% %   lam* = 13.919824 (60 nodes in s2); critical mode n = 1.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Road_Field_Berestycki
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
evalin('base','stability_form = ''P'';');  % form of the stability test

% Default parameters (lam = 0.9*lam*), overridden by the strings in params
Dr = 2;     d = 1;      mu = 1;     nu = 1;
lam = 0.9*13.919824;
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
r = pde_var(s1,[0,1]);
f = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(r,'t')==Dr*diff(r,s1,2)+nu*subs(f,s2,0)-mu*r;
         diff(f,'t')==d*diff(f,s1,2)+d*diff(f,s2,2)+lam*f;
         subs(r,s1,0)==0;           subs(r,s1,1)==0;
         subs(f,s1,0)==0;           subs(f,s1,1)==0;
         subs(diff(f,s2),s2,0)==(nu/d)*subs(f,s2,0)-(mu/d)*r;
         subs(f,s2,1)==0];


% Example-specific settings: the light tier, with the 1-D storage from the
% multiplier only (s-lengths 0 of both integrals) and the 2-D storage with
% theta-lengths [1 1] (multiplier) and [1 2] (integral)
st = lpisettings_sop('light');
st.lpi.pos.len_inc{1} = 0;
st.lpi.pos.len_com{1} = 0;
st.lpi.pos.len_sos{2} = [1 1;1 2;0 0];

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{berestycki2013influence,
%   title={The influence of a line with fast diffusion on Fisher-KPP propagation},
%   author={Berestycki, Henri and Roquejoffre, Jean-Michel and Rossi, Luca},
%   journal={Journal of Mathematical Biology},
%   volume={66},
%   pages={743--766},
%   year={2013},
%   note={arXiv:1210.3721; journal DOI 10.1007/s00285-012-0508-1 unverified}
% }
