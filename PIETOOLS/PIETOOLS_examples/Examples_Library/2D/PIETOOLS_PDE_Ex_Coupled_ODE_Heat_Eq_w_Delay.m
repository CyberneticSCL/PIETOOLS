function [PDE_t,st] = PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay(GUI,params) % MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_COUPLED_ODE_HEAT_EQ_W_DELAY(GUI,PARAMS)
% returns the delayed ODE-heat equation cascade below, the delays written
% as transport states x3(s2), x4(s1,s2), as a 'pde_struct' object, and
% the recommended settings ST of its stability test. It sets stability = 1
% and stability_form = 'Q' (generalized storage T'Q) in the base
% workspace; a nonzero GUI prints that no GUI representation exists.
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/21/2026
%               for the tests this file names                               % MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % ODE coupled to heat equation with delay (Kang, 2017):
% % ODE         X_{t} = A*X(t) + A1*X(t-tau) + B*x(t,s1=0)
% % PDE         x_{t} = x_{s1s1} + a*x + a2*x(t-tau);     s1 in [0,1]
% % With BCs    x_{s1}(t,s1=0) = 0;   
% %             x(t,s1=1) = 0;
% %
% % Let x1(t)=X(t) and x2(t,s1)=x(t,s1). Let further x3(t,s2) and x4(t,s1,s2)
% % be such that:
% %     x3(t,s2=1+tau)    = x1(t)     with x3_{t} = x3_{s2}   s2 in [1,1+tau]
% %     x4(t,s1,s2=1+tau) = x2(t,s1)  with x4_{t} = x4_{s2}
% % Then we can write the system as
% %
% % ODE:    x1_{t} = A*x1(t) + A1*x3(t,1) + B*x2(t,s1=0);
% % PDE:    x2_{t} = x2_{s1s1} + a*x2 + a2*x4(t,s1,s2=1);
% %         x3_{t} = x3_{s2};
% %         x4_{t} = x4_{s2};
% % BCs:    x2_{s1}(t,s1=0) = 0;
% %         x2(t,s1=1) = 0;
% %         x3(t,s2=1+tau) = x1(t);
% %         x4(t,s1,s2=1+tau) = x2(t,s1);
% %         x4_{s1}(t,s1=0,s2) = 0;
% %         x4(t,s1=1,s2) = 0;
% %
% % Parameters tau, A, A1, B, a and a2 can be set.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay
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
% Initial coding DJ - 08/22/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the
%   recommended settings of the stability test; stability_form = 'Q' set in
%   the base workspace; help revised (OUTPUT lists no PDE_b, which is not
%   returned), license added.
% function [PDE_t] = PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay(GUI,params) % MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/21/2026
%evalin('base','stability_dual = 1;');

% Default parameters, overridden by the assignment strings in params
ne = 1; tau = 1; A = -1; A1 = 0.5; B = -1; a = 1; a2 = 0.5;
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
pde_var x1 x2 x3 x4
x2.vars = s1;       x2.dom = [0,1];
x3.vars = s2;       x3.dom = [1,1+tau];
x4.vars = [s1;s2];  x4.dom = [0,1; 1,1+tau];
PDE_t = [diff(x1,'t')==A*x1 +B*subs(x2,s1,0)+A1*subs(x3,s2,1);
         diff(x2,'t')==diff(x2,s1,2)+a*x2+a2*subs(x4,s2,1);
         diff(x3,'t')==diff(x3,s2);
         diff(x4,'t')==diff(x4,s2);
         subs(diff(x2,s1),s1,0)==0;     subs(x2,s1,1)==0;
         subs(x3,s2,1+tau)==x1;
         subs(diff(x4,s1),s1,0)==0;     subs(x4,s1,1)==0;
         subs(x4,s2,1+tau)==x2];


% %%% Terms format
% % Initialize the state variables and input.
% PDE_t.x{1}.vars = [];
% PDE_t.x{2}.vars = s1;       PDE_t.x{2}.dom = [0,1];
% PDE_t.x{3}.vars = s2;       PDE_t.x{3}.dom = [1,1+tau];
% PDE_t.x{4}.vars = [s1;s2];  PDE_t.x{4}.dom = [0,1; 1,1+tau];
% 
% % ODE: x1_{t} = A*x1(t)
% PDE_t.x{1}.term{1}.x = 1;
% PDE_t.x{1}.term{1}.C = A;
% % ODE: x1_{t} = ... + A1*x3(t,1)
% PDE_t.x{1}.term{2}.x = 3;
% PDE_t.x{1}.term{2}.loc = 1;
% PDE_t.x{1}.term{2}.C = A1;
% % ODE: x1_{t} = ... + B*x2(t,s1=0);
% PDE_t.x{1}.term{3}.x = 2;
% PDE_t.x{1}.term{3}.loc = 0;
% PDE_t.x{1}.term{3}.C = B;
% 
% % PDE: x2_{t} = x2_{s1s1} + a*x2
% PDE_t.x{2}.term{1}.x = [2; 2];
% PDE_t.x{2}.term{1}.D = [2; 0];
% PDE_t.x{2}.term{1}.C = [1, a];
% % PDE: x2_{t} = ... + a2*x4(t,s1,s2=1);
% PDE_t.x{2}.term{2}.x = 4;
% PDE_t.x{2}.term{2}.loc = [s1,1];
% PDE_t.x{2}.term{2}.C = a2;
% 
% % PDE: x3_{t} = x3_{s2};
% PDE_t.x{3}.term{1}.x = 3;
% PDE_t.x{3}.term{1}.D = 1;
% 
% % PDE: x4_{t} = x4_{s2};
% PDE_t.x{4}.term{1}.x = 4;
% PDE_t.x{4}.term{1}.D = [0,1];
% 
% 
% % BC1: 0 = x2_{s1}(t,s1=0);
% PDE_t.BC{1}.term{1}.x = 2;
% PDE_t.BC{1}.term{1}.D = 1;
% PDE_t.BC{1}.term{1}.loc = 0;
% 
% % BC2: 0 = x2(t,s1=1)
% PDE_t.BC{2}.term{1}.x = 2;
% PDE_t.BC{2}.term{1}.loc = 1;
% PDE_t.BC{2}.term{1}.D = 0;
% 
% % BC3: 0 = x3(t,s2=1+tau) - x1(t);
% PDE_t.BC{3}.term{1}.x = 3;          PDE_t.BC{3}.term{2}.x = 1;
% PDE_t.BC{3}.term{1}.loc = 1+tau;    PDE_t.BC{3}.term{2}.C = -1;
% 
% % BC4: 0 = x4(t,s1,s2=1+tau) - x2(t,s1);
% PDE_t.BC{4}.term{1}.x = 4;              PDE_t.BC{4}.term{2}.x = 2;
% PDE_t.BC{4}.term{1}.loc = [s1,1+tau];   PDE_t.BC{4}.term{2}.loc = s1;          
%                                         PDE_t.BC{4}.term{2}.C = -1; 
% 
% % BC5: 0 = x4_{s1}(t,s1=0,s2);
% PDE_t.BC{5}.term{1}.x = 4;
% PDE_t.BC{5}.term{1}.D = [1,0];
% PDE_t.BC{5}.term{1}.loc = [0,s2];
% 
% % BC6: 0 = x4(t,s1=1,s2)
% PDE_t.BC{6}.term{1}.x = 4;          
% PDE_t.BC{6}.term{1}.loc = [1,s2];
% PDE_t.BC{6}.term{1}.D = [0,0];              
                                        
                                        
% Example-specific settings                                                 % MMP, 11/21/2026
% (stripped tier; slack lift degree at its minimum Dmin; Q of degree 0,
% except degree 1 in the output variable of its integral kernels)
st = lpisettings_sop('stripped');                                           % MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/21/2026
st.sop.deg_Q{2} = [0 0 1];                                                  % MMP, 11/21/2026

% No GUI representation of this system
if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{kang2017boundary,
%   title={Boundary control of delayed ODE--heat cascade under actuator saturation},
%   author={Kang, Wen and Fridman, Emilia},
%   journal={Automatica},
%   volume={83},
%   pages={252--261},
%   year={2017},
%   publisher={Elsevier}
% }