function [PDE_t, st] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params) % MMP, 11/18/2026
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
% - st:         (optional) example-specific settings structure (header)     % MMP, 11/18/2026
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
%
% %---------------------------------------------------------------------% %
% % 2D advection-diffusion equation (Holmes, 1994):
% % PDE         x_{t}   = C*(x_{s1s1} + x_{s2s2}) + b1*x_{s1} + b2*x_{s2}
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% %             x(s1=1) = 0;   x(s2=1) = 0;
% %
% % Parameters C, b1, and b2 can be set.
% %                                                                         % MMP, 11/18/2026
% % Example-specific settings (output st), for                              % MMP, 11/18/2026
% % PIETOOLS_stability_sop(PIE,st,'Q'): lpisettings_sop('light') with       % MMP, 11/18/2026
% % the slack lift degree at its minimum, st.lpi.ineq.dD = 0, and Q of      % MMP, 11/18/2026
% % degree 0 in two variables, st.sop.deg_Q{2} = 0. MOSEK, 11/18/2026,      % MMP, 11/18/2026
% % library defaults: 8,009 decision variables, 815 rows, OPT_STATUS        % MMP, 11/18/2026
% % 0.984, eta 1.6e-11; light tier: 72,100 variables, 3,281 rows            % MMP, 11/18/2026
% % (OPT_STATUS 1.01, eta 7.9e-10). Verified at the default parameters      % MMP, 11/18/2026
% % only.                                                                   % MMP, 11/18/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/18/2026: Second output st, example-specific settings that
%                   certify the library default at fewer decision
%                   variables than the light tier (header);
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq(GUI,params)   % MMP, 11/18/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

C = 1;  b1 = 0.5;  b2 = 2;  ne = 1;
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
x = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(x,'t')==C*(diff(x,s1,2)+diff(x,s2,2)) +b1*diff(x,s1)+b2*diff(x,s2);
         subs(x,s1,0)==0;   subs(x,s1,1)==0;
         subs(x,s2,0)==0;   subs(x,s2,1)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
% 
% % PDE: x_{t} = [C, C] * [x_{s1s1}; x_{s2s2}];
% PDE_t.x{1}.term{1}.D = [2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [C*eye(ne), C*eye(ne)];
% 
% % PDE: x_{t} = ... + [b1, b2] * [x_{s1}; x_{s2}];
% PDE_t.x{1}.term{2}.D = [1,0; 0,1];
% PDE_t.x{1}.term{2}.C = [b1*eye(ne), b2*eye(ne)];
% 
% % BC1: 0 = x(s1,0)                     % BC3: 0 = x(s1,1)
% PDE_t.BC{1}.term{1}.loc = [s1,0];      PDE_t.BC{3}.term{1}.loc = [s1,1];
% % BC2: 0 = x(0,s2)                     % BC4: 0 = x(1,s2)
% PDE_t.BC{2}.term{1}.loc = [0,s2];      PDE_t.BC{4}.term{1}.loc = [1,s2];


% Example-specific settings (header): light, slack lift degree at its       % MMP, 11/18/2026
% minimum (dD = 0), Q of degree 0 in two variables (st.sop.deg_Q)           % MMP, 11/18/2026
if nargout>1                                                                % MMP, 11/18/2026
    st = lpisettings_sop('light');                                          % MMP, 11/18/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
    st.sop.deg_Q{2} = 0;                                                    % MMP, 11/18/2026
end                                                                         % MMP, 11/18/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.