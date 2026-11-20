function [PDE_t, st] = PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang(GUI,params) % MMP, 11/17/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'b=1;','d=0.1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) settings structure that certifies the example    % MMP, 11/17/2026
%               (application-specific settings, settings report             % MMP, 11/17/2026
%               11/16/2026 Sec. 5.8).                                       % MMP, 11/17/2026
%
% %---------------------------------------------------------------------% %
% % Age-structured population with spatial diffusion (linear McKendrick
% % model with diffusion; Kang 2022) [standard linear form, unverified],
% % s1 = age a, s2 = space x:
% % PDE         q_{t}(t,s1,s2) = -q_{s1} + d*q_{s2s2} - mu*q;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    q(t,s1=0,s2) = b*int_0^1 q(t,s1,s2) ds1;   (renewal)
% %             q(t,s1,s2=0) = 0;   q(t,s1,s2=1) = 0;
% % q is first order in s1 and second order in s2.
% %
% % Parameters b, d and mu can be set; default d = 0.1, mu = 1,
% % b = 0.9*b*.
% %
% % Threshold [derived from the modes, exact]: in the sine modes in s2 the
% % first mode decays iff R_eff = b*(1-exp(-c))/c < 1, c = mu + d*pi^2:
% %   stable iff b < b* = c/(1-exp(-c)) = 2.302685   (c = 1.986960).
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Second output st, the settings that certify the
%                   example in the storage form 'P', which the light
%                   tier does not (maintainer: application-specific
%                   settings where the tiers do not); header;
% function [PDE_t] = PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang(GUI,params) % MMP, 11/17/2026 (was)
% MMP, 11/18/2026: st for the form 'Q' of the library: light with the
%                   slack lift degree at its minimum, 12.6 times fewer
%                   variables than the 11/17 st in 'Q' (header);
% MMP, 11/19/2026: st updated to the per-example search of 11/19/2026
%                   (maintainer: no default settings): st.sop.deg_Q{2} = [0 1 1];
%                   the form of its tests named in the workspace.
% MMP, 11/19/2026: the settings notes of the header (11/17-11/19) removed: they
%                   recorded tests and files the user does not have, not the problem;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% MMP, 11/19/2026: the dated test records of the header removed: they
%                   recorded tests and files the user does not have, not
%                   the problem;

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026

d = 0.1;    mu = 1;     b = 0.9*2.302685;
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
q = pde_var([s1;s2],[0,1;0,1],[1;2]);           % 1st order in s1 (age), 2nd in s2
PDE_t = [diff(q,'t')==-diff(q,s1)+d*diff(q,s2,2)-mu*q;
         subs(q,s1,0)==b*int(q,s1,[0,1]);
         subs(q,s2,0)==0;           subs(q,s2,1)==0];


% Application-specific settings (header): the heavy tier, for the           % MMP, 11/18/2026 (was); MMP, 11/17/2026
% storage form 'P', without the slack weight raise in two variables         % MMP, 11/18/2026 (was); MMP, 11/17/2026
% Example-specific settings                                                 % MMP, 11/19/2026
if nargout>1                                                                % MMP, 11/17/2026
%    st = lpisettings_sop('heavy');                                         % MMP, 11/18/2026 (was); MMP, 11/17/2026
    st = lpisettings_sop('light');                                          % MMP, 11/18/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
%    st.sop.dw_feas(2) = 0;                                                 % MMP, 11/18/2026 (was); MMP, 11/17/2026
    st.sop.deg_Q{2} = [0 1 1];                                              % MMP, 11/19/2026
end                                                                         % MMP, 11/17/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{kang2022effects,
%   title={Effects of diffusion on the principal eigenvalue for age-structured models with random diffusion},
%   author={Kang, Hao},
%   journal={Proceedings of the Royal Society of Edinburgh Section A},
%   volume={152},
%   number={1},
%   pages={258--280},
%   year={2022},
%   doi={10.1017/prm.2021.10}
% }
