function [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'r1=10;','r2=15;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Cascade of a 1D and a 2D heat equation with variable coefficients:
% % the 1D state x1 is the Dirichlet datum of the 2D state x2 on the edge
% % s2 = 0.
% % PDE         x1_{t}(t,s1)    = (1+s1)*x1_{s1s1} + r1*x1;   s1 in [0,1]
% %             x2_{t}(t,s1,s2) = x2_{s1s1} + (1+s2^2)*x2_{s2s2} + r2*x2;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    x1(t,s1=0) = 0;   x1(t,s1=1) = 0;
% %             x2 = 0 at s1 = 0, 1;
% %             x2(t,s1,s2=0) = x1(t,s1);   x2(t,s1,s2=1) = 0;
% %
% % Parameters r1 and r2 can be set; default r1 = 0.9*r1*, r2 = 0.9*r2*.
% %
% % Thresholds of the blocks [finite differences, second-order central,
% % N interior nodes per direction]:
% %   r1* = 14.504004 (N = 40), 14.507976 (N = 60);
% %   r2* = 22.199946 (N = 40), 22.205980 (N = 60).
% % The cascade is triangular, hence stable iff r1 < r1* and r2 < r2*
% % [inferred, not computed for the coupled system].
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier certifies the example, so
%                   no second output st (application-specific settings);
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
%                   the entry of 11/17/2026 (no second output st) no longer
%                   holds;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% MMP, 11/19/2026: the dated test records of the header removed: they
%                   recorded tests and files the user does not have, not
%                   the problem;
% function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff(GUI,params)  % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026

r1 = 0.9*14.507976;     r2 = 0.9*22.205980;
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
x1 = pde_var(1,s1,[0,1]);
x2 = pde_var(1,[s1;s2],[0,1;0,1]);
PDE_t = [diff(x1,'t')==(1+s1)*diff(x1,s1,2)+r1*x1;
         diff(x2,'t')==diff(x2,s1,2)+(1+s2^2)*diff(x2,s2,2)+r2*x2;
         subs(x1,s1,0)==0;          subs(x1,s1,1)==0;
         subs(x2,s1,0)==0;          subs(x2,s1,1)==0;
         subs(x2,s2,0)==x1;         subs(x2,s2,1)==0];


% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('light');                                              % MMP, 11/19/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/19/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
