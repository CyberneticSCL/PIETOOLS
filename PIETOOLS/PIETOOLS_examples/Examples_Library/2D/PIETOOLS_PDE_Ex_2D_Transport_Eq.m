function [PDE_t, st] = PIETOOLS_PDE_Ex_2D_Transport_Eq(GUI,params)          % MMP, 11/18/2026
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
% % 2D transport equation
% % PDE         x_{t}  = c1*x_{s1} + c2*x_{s2}
% % With BCs    x(s1=0) = 0;   x(s2=0) = 0;
% %
% % Parameters c1 and c2 can be set.
% %---------------------------------------------------------------------% %
%
% MMP, 11/18/2026: Second output st, example-specific settings that
%                   certify the library default at fewer decision
%                   variables than the light tier (header);
% MMP, 11/19/2026: st updated to the per-example search of 11/19/2026
%                   (maintainer: no default settings): st.lpi.ineq.psatz{2} = 'none';
%                   the form of its tests named in the workspace.
% MMP, 11/19/2026: the settings notes of the header removed: they recorded
%                   tests and files the user does not have, not the problem;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_Transport_Eq(GUI,params)            % MMP, 11/18/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026
%evalin('base','stability_dual = 1;');

c1 = -1;     c2 = -1;     ne = 1;
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
PDE_t = [diff(x,'t')==c1*diff(x,s1)+c2*diff(x,s2);
         subs(x,s1,0)==0;       subs(x,s2,0)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
% 
% % PDE: x_{t} = [c1, c2] * [x_{s1}; x_{s2}]
% PDE_t.x{1}.term{1}.D = [1,0; 0,1];
% PDE_t.x{1}.term{1}.C = [c1*eye(ne), c2*eye(ne)];
% 
% % BC1: 0 = x(s1,0)
% PDE_t.BC{1}.term{1}.loc = [s1,0];
% % BC2: 0 = x(0,s2)
% PDE_t.BC{2}.term{1}.loc = [0,s2];


% Example-specific settings                                                 % MMP, 11/19/2026
if nargout>1                                                                % MMP, 11/18/2026
    st = lpisettings_sop('light');                                          % MMP, 11/18/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
    st.sop.deg_Q{2} = 0;                                                    % MMP, 11/18/2026
    st.lpi.ineq.psatz{2} = 'none';                                          % MMP, 11/19/2026
end                                                                         % MMP, 11/18/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end