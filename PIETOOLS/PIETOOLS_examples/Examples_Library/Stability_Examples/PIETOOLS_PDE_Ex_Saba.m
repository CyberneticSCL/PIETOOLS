function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Saba(index,GUI,params)          % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_SABA(INDEX,GUI,PARAMS) returns the
% 2-state linear hyperbolic system below with the parameters selected by
% INDEX, in the term-based and in the batch format, and the recommended
% settings ST of its stability test. It sets stability = 1 and
% stability_form = 'P' (storage <T x, P T x>, P >= eppos I) in the base
% workspace. With GUI = 1 it also loads the associated GUI save file into
% a new PIETOOLS_PDE_GUI window.
%
% INPUT
% - index:      4 or 4.1, 4.2, 4.3, 4.4, 4.5: the parameter sets of the
%               body; any other value raises an error.
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
%               for the test this file names, by case                       % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Saba 2019 examples (see reference below): 
% % PDE        x1_{t} = sig1*x2 - (1/r1)*x1_{s} 
% %            x2_{t} = sig2*x1 + (1/r2)*x2_{s}
% % with BCs   x1(s=0) = qb*x2(s=0)
% %            x2(s=1) = pb*x1(s=1)
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Saba
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
% MMP, 11/19/2026: Output st, the example-specific settings of each case,
%                   and the form of its test named in the workspace
%                   (stability_form), from the per-example search of
%                   11/19/2026 (maintainer: no default settings);
% MMP, 11/19/2026: case 4.5 takes sig2 = 1.048, the value of its source (Peet,
%                   Automatica 2021, Sec. 10, Ex. 5: stability for sig2 <=
%                   1.048 at pb = -0.4); 2.048, a typo, was unstable (real
%                   eigenvalue 1.682); maintainer: "Apply all four";
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Saba(index,GUI,params)           % MMP, 11/19/2026 (was)
% MMP, 11/19/2026: Release comment pass (comment-only).

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;');
% evalin('base','stability_dual = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/19/2026


%%%% There are several examples included here. Add a decimal to
%%%% your example input to specify a particular set of parameters.
if index==4 || index==4.1
 r1=.8;   r2=1.1;   sig1=2.3;   sig2=-3.5;    qb=-.7;   pb=.5;
 dir = fullfile(root,'PIETOOLS_PDE_Ex_Saba_1_GUI.mat');
 
elseif index==4.2
 r1=.5;   r2=1.1;   sig1=1;     sig2=-.1;     qb=1.2;   pb=0;
 dir = fullfile(root,'PIETOOLS_PDE_Ex_Saba_2_GUI.mat');
 
elseif index==4.3
 r1=.5;   r2=1.1;   sig1=1;     sig2=1;       qb=1.2;   pb=-.7;
 dir = fullfile(root,'PIETOOLS_PDE_Ex_Saba_3_GUI.mat');
 
elseif index==4.4
 r1=.5;   r2=1.1;   sig1=1;     sig2=.663;    qb=1.2;   pb=0;   % max sig2=.663 
 dir = fullfile(root,'PIETOOLS_PDE_Ex_Saba_4_GUI.mat');
 
elseif index==4.5
% r1=.5;   r2=1.1;   sig1=1;     sig2=2.048;   qb=1.2;   pb=-.4; % max sig2=1.049 % MMP, 11/19/2026 (was)
 r1=.5;   r2=1.1;   sig1=1;     sig2=1.048;   qb=1.2;   pb=-.4; % max sig2=1.049 % MMP, 11/19/2026
 dir = fullfile(root,'PIETOOLS_PDE_Ex_Saba_5_GUI.mat');
 
else
 error(['Example ',num2str(indx),'does not exist, please choose from 4.1, 4.2, 4.3, 4.4 or 4.5']);
end

npars = length(params);
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end


 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [0 sig1; sig2 0]; 
 PDE_b.A1 = [-1/r1 0; 0 1/r2];
 PDE_b.B = [1 -qb, 0, 0; 0, 0, -pb, 1];

 
 %%% Term-based input format
 PDE_t.x{1}.vars = s;
 PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = [0,sig1;sig2,0] * x
 PDE_t.x{1}.term{1}.C = [0, sig1; sig2, 0];

 % PDE: x_{t} = ... + [-1/r1,0;0,1/r2] * x_{s}
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = [-1/r1, 0; 0, 1/r2];
 
 % BCs: 0 = [1,-qb;0,0] * x(0)
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [1 -qb; 0 0];

 % BCs: 0 = ... + [0,0;-pb,1] * x(1)
 PDE_t.BC{1}.term{2}.loc = 1;
 PDE_t.BC{1}.term{2}.C = [0 0; -pb 1];

% Example-specific settings, by case                                        % MMP, 11/19/2026
% (stripped tier; 4.1: 1-D storage s-lengths 1 (lower and upper integral);
% 4.2: multiplier only, theta-length 1, slack lift degree at its minimum
% Dmin; 4.3, 4.5: slack lift degree Dmin; 4.4: upper-integral s-length 1,
% theta-lengths 2, 1, 1 (multiplier, upper, lower integral), Dmin)
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
if index==4 || index==4.1                                                   % MMP, 11/19/2026
    st.lpi.pos.len_inc{1} = 1;                                              % MMP, 11/19/2026
    st.lpi.pos.len_com{1} = 1;                                              % MMP, 11/19/2026
elseif index==4.2                                                           % MMP, 11/19/2026
    st.lpi.pos.len_inc{1} = 0;                                              % MMP, 11/19/2026
    st.lpi.pos.len_com{1} = 0;                                              % MMP, 11/19/2026
    st.lpi.pos.len_sos{1} = [1;2;2];                                        % MMP, 11/19/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/19/2026
elseif index==4.3 || index==4.5                                             % MMP, 11/19/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/19/2026
elseif index==4.4                                                           % MMP, 11/19/2026
    st.lpi.pos.len_com{1} = 1;                                              % MMP, 11/19/2026
    st.lpi.pos.len_sos{1} = [2;1;1];                                        % MMP, 11/19/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/19/2026
end                                                                         % MMP, 11/19/2026
 
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(dir);
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{saba2019stability,
%   title={Stability Analysis for a Class of Linear 2x2 Hyperbolic PDEs Using a Backstepping Transform},
%   author={Saba, David Bou and Argomedo, Federico Bribiesca and Auriol, Jean and Di Loreto, Michael and Di Meglio, Florent},
%   journal={IEEE Transactions on Automatic Control},
%   year={2019},
%   publisher={IEEE}
% }