function [PDE_t,st] = PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'L=0.5;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Elastic-thermoelastic transmission problem, 1D (Marzocchi, Munoz
% % Rivera, Naso 2002) [unverified standard form], a = b = k = alpha =
% % beta = 1:
% % PDE         u_{tt} = u_{xx}                         x in [0,L]
% %             v_{tt} = v_{xx} - theta_{x},
% %             theta_{t} = theta_{xx} - v_{xt}         x in [L,2L]
% % Transmission at x = L: continuity of displacement and stress.
% %
% % Let e1 = u_{x}, p = u_{t} on [0,L], and e2 = v_{x}, q = v_{t}, th =
% % theta on [L,2L], shifted to s1 in [0,L]:
% % PDE         e1_{t} = p_{s1};   p_{t} = e1_{s1};
% %             e2_{t} = q_{s1};   q_{t} = e2_{s1} - th_{s1};
% %             th_{t} = th_{s1s1} - q_{s1};
% % With BCs    e1(t,s1=0) = 0;   q(t,s1=L) = 0;   th(t,0) = th(t,L) = 0;
% %             p(t,s1=L) = q(t,s1=0);   e1(t,s1=L) = e2(t,s1=0);
% % Modification: a free end at x = 0 (e1(0) = 0) replaces the clamped end
% % of the paper, which leaves the uniform-strain mode as an equilibrium of
% % the Riemann form (a zero eigenvalue of the PIE).
% %
% % Parameter L can be set; default L = 0.5.
% %
% % The paper shows exponential decay of the energy for the clamped
% % problem [secondary citations and abstract]. Exponential stability of
% % the modified system is not established here; no threshold.
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
% function [PDE_t] = PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi(GUI,params)  % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''P'';');  % form of the stability test     % MMP, 11/19/2026

L = 0.5;
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
e1 = pde_var(s1,[0,L],1);   pp = pde_var(s1,[0,L],1);
e2 = pde_var(s1,[0,L],1);   q = pde_var(s1,[0,L],1);
th = pde_var(s1,[0,L]);
PDE_t = [diff(e1,'t')==diff(pp,s1);
         diff(pp,'t')==diff(e1,s1);
         diff(e2,'t')==diff(q,s1);
         diff(q,'t')==diff(e2,s1)-diff(th,s1);
         diff(th,'t')==diff(th,s1,2)-diff(q,s1);
         subs(e1,s1,0)==0;          subs(q,s1,L)==0;
         subs(th,s1,0)==0;          subs(th,s1,L)==0;
         subs(pp,s1,L)==subs(q,s1,0);
         subs(e1,s1,L)==subs(e2,s1,0)];


% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.lpi.pos.len_com{1} = 0;                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [0;2;1];                                            % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{marzocchi2002asymptotic,
%   title={Asymptotic behaviour and exponential stability for a transmission problem in thermoelasticity},
%   author={Marzocchi, Alfredo and Mu{\~n}oz Rivera, Jaime E and Naso, Maria Grazia},
%   journal={Mathematical Methods in the Applied Sciences},
%   volume={25},
%   number={11},
%   pages={955--980},
%   year={2002},
%   doi={10.1002/mma.323}
% }
