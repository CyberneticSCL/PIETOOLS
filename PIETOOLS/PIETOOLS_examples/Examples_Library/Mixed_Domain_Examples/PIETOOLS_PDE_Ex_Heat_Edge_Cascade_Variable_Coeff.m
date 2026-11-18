function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff(GUI,params)
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
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, MOSEK, 11/16/2026), at
% % (r1,r2) = f*(r1*,r2*), N = 60 values:
% %   - f = 0.9: UNKNOWN with the stock light settings (508,228 decision
% %     variables, after one raise); certified with the proposed light row
% %     of the settings proposal of 11/15/2026 (510,057) and with the
% %     stock heavy settings (1,023,444).
% %   - bisection of f on [0.5,1], 8 steps, light: stock certifies
% %     f = 0.73633, the proposed light row f = 0.97656.
% %                                                                         % MMP, 11/17/2026
% % No second output st: the storage of the light tier is the proposed      % MMP, 11/17/2026
% % light row since 11/17/2026, which certifies this example at f = 0.9     % MMP, 11/17/2026
% % (record above); the slack of 11/17/2026 (wcross 0, 2-D faces at         % MMP, 11/17/2026
% % offset 1) was not rerun.                                                % MMP, 11/17/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier certifies the example, so
%                   no second output st (application-specific settings);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');

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


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
