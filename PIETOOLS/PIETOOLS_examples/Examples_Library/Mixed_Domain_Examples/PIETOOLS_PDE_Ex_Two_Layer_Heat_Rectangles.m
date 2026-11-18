function [PDE_t] = PIETOOLS_PDE_Ex_Two_Layer_Heat_Rectangles(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=10;','H=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % Two-layer heat conduction with imperfect contact (Hickson, Barry,
% % Mercer 2009, 1D problem) on two rectangles [derived]: layer 1 on
% % [-1,0]x[0,1], layer 2 on [0,1]x[0,1], interface x = 0. Layer 1 is
% % reflected (x -> -s1) onto the unit square, so the interface is the
% % edge s1 = 0 of both states:
% % PDE         u1_{t} = k1*(u1_{s1s1} + u1_{s2s2}) + lam*u1;
% %             u2_{t} = k2*(u2_{s1s1} + u2_{s2s2}) + lam*u2;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    -k1*u1_{s1}(t,0,s2) = H*(u2(t,0,s2) - u1(t,0,s2));
% %              k2*u2_{s1}(t,0,s2) = H*(u2(t,0,s2) - u1(t,0,s2));
% %             u1 = u2 = 0 at s1 = 1, s2 = 0, s2 = 1.
% %
% % Parameters lam, k1, k2 and H can be set; default k1 = 1, k2 = 2,
% % H = 1, lam = 0.9*lam*.
% %
% % Threshold [derived, finite differences]: sine modes in s2 (mode n adds
% % -k_i*(n*pi)^2), finite differences of the 1D transmission problem
% % with ghost-node interface conditions, modes n = 1..3, default k1, k2, H:
% %   lam* = 13.686807 (N = 40), 13.686938 (N = 60).
% % The paper treats the 1D problem and states no stability threshold.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), lam = 0.9*lam*: certified after one raise with
% % the stock light settings (751,936 decision variables) and with the
% % proposed light row of the settings proposal of 11/15/2026 (754,180).
% %                                                                         % MMP, 11/17/2026
% % No second output st: the light tier certifies this example              % MMP, 11/17/2026
% % (record above). Its storage is the proposed light row since             % MMP, 11/17/2026
% % 11/17/2026; its slack of 11/17/2026 (wcross 0, 2-D faces at offset      % MMP, 11/17/2026
% % 1) was not rerun.                                                       % MMP, 11/17/2026
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

k1 = 1;     k2 = 2;     H = 1;
lam = 0.9*13.686938;
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
u1 = pde_var([s1;s2],[0,1;0,1]);
u2 = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(u1,'t')==k1*diff(u1,s1,2)+k1*diff(u1,s2,2)+lam*u1;
         diff(u2,'t')==k2*diff(u2,s1,2)+k2*diff(u2,s2,2)+lam*u2;
         subs(diff(u1,s1),s1,0)==-(H/k1)*(subs(u2,s1,0)-subs(u1,s1,0));
         subs(diff(u2,s1),s1,0)==(H/k2)*(subs(u2,s1,0)-subs(u1,s1,0));
         subs(u1,s1,1)==0;          subs(u2,s1,1)==0;
         subs([u1;u2],s2,0)==0;     subs([u1;u2],s2,1)==0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{hickson2009critical,
%   title={Critical times in multilayer diffusion. Part 1: Exact solutions},
%   author={Hickson, Roslyn I and Barry, Steven I and Mercer, Geoffry N},
%   journal={International Journal of Heat and Mass Transfer},
%   volume={52},
%   pages={5776--5783},
%   year={2009},
%   doi={10.1016/j.ijheatmasstransfer.2009.08.013}
% }
% (Volume and pages unverified.)
