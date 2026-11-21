function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(GUI,params) % DB, MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_PARABOLIC_EQ_GAHLAWAT(GUI,PARAMS)
% returns the parabolic equation with spatially varying coefficients below
% (Gahlawat 2017, Problem 2), in the term-based and in the batch format,
% and the recommended settings ST of its stability tests. It sets
% stability = 1, stability_dual = 1 and stability_form = 'Q' (generalized
% storage T'Q; dual form 'Qd') in the base workspace. With GUI = 1 it also
% loads the associated GUI save file into a new PIETOOLS_PDE_GUI window.
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
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % DB, MMP, 11/21/2026
%               for the tests this file names                               % DB, MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Diffusion Equation Problem 2 from Gahlawat 2017:
% % PDE        x_{t} = c(s)*x + b(s)*x_{s} + a(s)*x_{ss} 
% % with BCs   x(s=0) = 0
% %            x_{s}(s=1) = 0
% %
% % and with   a(s) = s^3 - s^2 + 2
% %            b(s) = 3*s^2 - 2*s
% %            c(s) = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam
% %
% % Parameters lam, ne (state size) and the coefficients a, b and c can be
% % set; c is computed from lam before PARAMS is evaluated.
% % Unstable for lam > 4.66.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat
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
% DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: default lam = 4.65
%   (was 4.6); new output st, the recommended settings of the stability
%   tests; stability_form = 'Q' set in the base workspace; help and
%   comments revised.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat(GUI,params) % DB, MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')
evalin('base','stability_form = ''Q'';');  % form of the stability test     % DB, MMP, 11/21/2026

% Specify the parameters
ne = 1; lam = 4.65;
a = s^3 - s^2 + 2;
b = 3*s^2 - 2*s;
c = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;
 
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = ne;
PDE_b.dom = [0,1];

PDE_b.A2 = a*eye(ne);
PDE_b.A1 = b*eye(ne);
PDE_b.A0 = c*eye(ne);

PDE_b.B = [eye(ne), zeros(ne), zeros(ne), zeros(ne);
    zeros(ne), zeros(ne), zeros(ne), eye(ne)];


%%% Term-based input format
PDE_t.x{1}.vars = s;   PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = cfun * x
PDE_t.x{1}.term{1}.C = c*eye(ne);

% PDE: x_{t} = ... + bfun * x_{s}
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.C = b*eye(ne);

% PDE: x_{t} = ... + afun * x_{ss}
PDE_t.x{1}.term{3}.D = 2;
PDE_t.x{1}.term{3}.C = a;

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;

% BCs: 0 = x_{s}(1)
PDE_t.BC{2}.term{1}.D = 1;
PDE_t.BC{2}.term{1}.loc = 1;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Parabolic_Eq_Gahlawat_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % DB, MMP, 11/21/2026
% (light tier; slack integral bases capped at total degree theta+s <=
% Dmax-1)
st = lpisettings_sop('light');                                              % DB, MMP, 11/21/2026
st.lpi.ineq.cap_inc = 'auto';                                               % DB, MMP, 11/21/2026

end
% @article{gahlawat2016convex,
%   title={A convex sum-of-squares approach to analysis, state feedback and output feedback control of parabolic PDEs},
%   author={Gahlawat, Aditya and Peet, Matthew M},
%   journal={IEEE Transactions on Automatic Control},
%   volume={62},
%   number={4},
%   pages={1636--1651},
%   year={2016},
%   publisher={IEEE}
% }