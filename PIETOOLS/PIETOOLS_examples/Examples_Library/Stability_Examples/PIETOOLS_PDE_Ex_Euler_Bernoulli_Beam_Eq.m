function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq(index,GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_EULER_BERNOULLI_BEAM_EQ(INDEX,GUI,
% PARAMS) returns the cantilevered Euler-Bernoulli beam below, in the
% state-space form INDEX selects, and the recommended settings ST of its
% stability tests. It sets stability = 1 and stability_dual = 1 in the
% base workspace, and stability_form = 'P' (storage <T x, P T x>,
% P >= eppos I) for INDEX 13, 13.1 or 'Q' (generalized storage T'Q) for
% INDEX 13.2; the dual tests take the forms 'Pd', 'Qd'. With GUI = 1
% (INDEX 13, 13.1) it also loads the associated GUI save file into a new
% PIETOOLS_PDE_GUI window.
%
% INPUT
% - index:      13 or 13.1: states x1 = u_{t}, x2 = u_{ss}, in the batch
%               and in the term-based format; 13.2: the state u, 4th
%               order in s, declared with pde_var (PDE_b then equals
%               PDE_t); any other value raises an error.
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names, by case                      % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Euler-Bernoulli beam equation from Peet (Example 8.1.0.1)
% % PDE        u_{tt} = -c*u_{ssss}
% % with BCs   u(s=0) = 0
% %            u_{s}(s=0) = 0 
% %            u_{ss}(s=1) = 0 
% %            u_{sss}(s=1) = 0 
% % 
% % Parameter c can be set.
% % For c > 0 the energy is conserved: stable, not exponentially stable.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq
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
% MMP, 11/19/2026: Output st, the example-specific settings of each case,
%                   and the form of its tests named in the workspace
%                   (stability_form), from the per-example search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq(index,GUI,params)  % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_dual = 1;')

% Specify the parameters
% Default parameter, overridden by the assignment strings in params
c = 0.1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

%%%% Two state-space representations are implemented. Add a decimal to
%%%% your example input to specify a particular set of parameters.
if index==13 || index==13.1
% % This case uses PDE states x1 = u_{t}, x2 = u_{ss}, so that the system becomes
% % PDE        x1_{t} = -c * x2_{ss}
% %            x2_{t} = x1_{ss}
% % with BCs   x1(0) = 0, x2(1) = 0, x1_{s}(0) = 0, x2_{s}(1) = 0
% Construct the PDE.
%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = 0;   PDE_b.n2 = 2;
PDE_b.dom = [0,1];

PDE_b.A2 = [0 -c; 1 0];

PDE_b.B = [ 1 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0;
    0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 1];


%%% Term-based input format
PDE_t.x{1}.vars = s;    PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = [0,-c;1,0] * x_{ss}
PDE_t.x{1}.term{1}.D = 2;
PDE_t.x{1}.term{1}.C = [0, -c; 1, 0];

% BCs: 0 = x1(0)
PDE_t.BC{1}.term{1}.loc = 0;
PDE_t.BC{1}.term{1}.C = [1,0];

% BCs: 0 = x2(1)
PDE_t.BC{2}.term{1}.loc = 1;
PDE_t.BC{2}.term{1}.C = [0,1];

% BCs: 0 = x1_{s}(0)
PDE_t.BC{3}.term{1}.D = 1;
PDE_t.BC{3}.term{1}.loc = 0;
PDE_t.BC{3}.term{1}.C = [1,0];

% BCs: 0 = x2_{s}(1)
PDE_t.BC{4}.term{1}.D = 1;
PDE_t.BC{4}.term{1}.loc = 1;
PDE_t.BC{4}.term{1}.C = [0,1];

%  % Note that BCs can also be coupled together:
%
%  PDE_t.BC{1}.term{1}.loc = 0;
%  PDE_t.BC{1}.term{1}.C = [1,0;0,0];
%
%  PDE_t.BC{1}.term{2}.loc = 1;
%  PDE_t.BC{1}.term{2}.C = [0,0;0,1];
%
%  PDE_t.BC{2}.term{1}.D = 1;
%  PDE_t.BC{2}.term{1}.loc = 0;
%  PDE_t.BC{2}.term{1}.C = [1,0;0,0];
%
%  PDE_t.BC{2}.term{2}.D = 1;
%  PDE_t.BC{2}.term{2}.loc = 1;
%  PDE_t.BC{2}.term{2}.C = [0,0;0,1];
% command-line formaty
% PDE_c=sys();
% x1 = state('pde');x2 = state('pde');
% eq_PDEs = [diff(x1,t) == -c*diff(x2,s,2); diff(x2,t)==diff(x1,s,2)];
% eq_BCs = [subs(x1,s,0) == 0; subs(x2,s,1) == 0; subs(diff(x1,s),s,0) ==0; subs(diff(x2,s),s,1)==0];
% PDE_c = addequation(PDE_c,[eq_PDEs;eq_BCs]);
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end
elseif index==13.2
     % In this case, the PDE states are set as u and u_{t}
    pvar s t; % define independent variables
    % Define dependent variables and system variable
    %   PDE: u_{tt} = -c*u_{ssss}                       | c = 0.1
    %   BCs: u(s=0) = 0,        u_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
    %        u_{s}(s=0) = 0,    u_{sss}(s=1) = 0        |  
    x=pde_var(s,[0,1]); 
    %% Define equations
    eq_PDE = diff(x,t,2)==-c*diff(x,s,4); % 	PDE: u_{tt} = -c*u_{ssss}
    eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,0)==0;
             subs(diff(x,s,2),s,1)==0;subs(diff(x,s,3),s,1)==0]; %   BCs: u(s=0) = 0, u_{ss}(s=1) = 0, u_{s}(s=0) = 0, u_{sss}(s=1) = 0
    %% initialize pde system;
    PDE_b = initialize([eq_PDE;eq_BC]);
    PDE_t = PDE_b;
    
    else
     error(['Example ',num2str(indx),'does not exist, please choose from 13.1 and 13.2']);
    end

% Example-specific settings and the form of the stability tests,            % MMP, 11/19/2026
% by case                                                                   % MMP, 11/19/2026
% (stripped tier with the slack lift degree at its minimum Dmin; for 13.2
% Q of degree 0 and the form 'Q'; otherwise the 1-D storage from the
% multiplier only, theta-length 1 (both integral s-lengths 0), and 'P')
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026
if index==13.2                                                              % MMP, 11/19/2026
    st.sop.deg_Q{1} = [0 0 0];                                              % MMP, 11/19/2026
    evalin('base','stability_form = ''Q'';');                               % MMP, 11/19/2026
else                                                                        % MMP, 11/19/2026
    st.lpi.pos.len_inc{1} = 0;                                              % MMP, 11/19/2026
    st.lpi.pos.len_com{1} = 0;                                              % MMP, 11/19/2026
    st.lpi.pos.len_sos{1} = [1;2;2];                                        % MMP, 11/19/2026
    evalin('base','stability_form = ''P'';');                               % MMP, 11/19/2026
end                                                                         % MMP, 11/19/2026

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