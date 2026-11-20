function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params)             % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % KISS Model from Holmes (1994):
% % PDE         x_{t}   = lam*x + c1*x_{s1s1} + c2*x_{s2s2}
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% %             x(s1=1) = 0;   x(s2=1) = 0;
% %
% % Parameters c1, c2 and lam can be set.
% % Stable for lam <= 2*pi^2 when c1=c2=1.
% %---------------------------------------------------------------------% %
%
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_KISS_Model(GUI,params)              % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026
%evalin('base','stability_dual = 1;');

c1 = 1; c2 = 1; lam = 19;   ne = 1;
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
x = pde_var(ne,[s1;s2],[0,1;0,1]);
PDE_t = [diff(x,'t')==lam*x+c1*diff(x,s1,2)+c2*diff(x,s2,2);
         subs(x,s1,0)==0;   subs(x,s1,1)==0;
         subs(x,s2,0)==0;   subs(x,s2,1)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
% 
% % PDE: x_{t} = [lam, c1, c2] * [x; x_{s1s1}; x_{s2s2}]
% PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [lam*eye(ne), c1*eye(ne), c2*eye(ne)];
% 
% % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
% PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
% PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];


% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.sop.deg_Q{2} = [0 1 1];                                                  % MMP, 11/19/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.