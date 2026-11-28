function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Transport_Eq(GUI,params)        % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_TRANSPORT_EQ(GUI,PARAMS) returns the
% transport equation below, in the term-based and in the batch format,
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
% - st:         example-specific settings (lpisettings structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Transport equation: 
% % PDE:        x_{t}  = v*x_{s}
% % with BC     x(s=0) = 0
% %
% % Parameter v (default -1) can be set. For v < 0 the BC is an inflow
% % condition and x(t) = 0 for t >= 1/|v|.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Transport_Eq
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
% MMP, 11/21/2026: Changes since PIETOOLS 2025: new output st, the
%   recommended settings of the stability tests; stability_form = 'Q' set
%   in the base workspace; help and comments revised.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Transport_Eq(GUI,params)         % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/21/2026
% Specify the parameters
v=-1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 1;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A1= v;

on = eye(PDE_b.n1);   ze = zeros(PDE_b.n1);
PDE_b.B = [on ze];


%%% Term-based input format
% Single 1D state component
PDE_t.x{1}.vars = s;
PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = x_{s}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C= v;
PDE_t.x{1}.term{1}.D = 1;

% BCs: 0 = x(0)
PDE_t.BC{1}.term{1}.loc = 0;


if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Transport_Eq_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % MMP, 11/21/2026
% (extreme tier; slack lift degree at its minimum Dmin; slack integral
% bases capped at total degree theta+s <= Dmax-1)
st = lpisettings('extreme');                                                % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/21/2026

end