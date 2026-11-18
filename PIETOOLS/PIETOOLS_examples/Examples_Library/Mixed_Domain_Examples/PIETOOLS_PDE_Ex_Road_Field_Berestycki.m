function [PDE_t] = PIETOOLS_PDE_Ex_Road_Field_Berestycki(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=10;','Dr=2;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % Road-field model (Berestycki, Roquejoffre, Rossi 2013, eq. (1.1),
% % linearized, on a box [derived]): a road r(s1) on the edge s2 = 0 of a
% % field f(s1,s2):
% % PDE         r_{t}(t,s1)    = Dr*r_{s1s1} + nu*f(t,s1,s2=0) - mu*r;
% %                                               s1 in [0,1]
% %             f_{t}(t,s1,s2) = d*(f_{s1s1} + f_{s2s2}) + lam*f;
% %                                               (s1,s2) in [0,1]^2
% % With BCs    r = 0 at s1 = 0, 1;   f = 0 at s1 = 0, 1 and at s2 = 1;
% %             -d*f_{s2}(t,s1,0) = mu*r(t,s1) - nu*f(t,s1,0);
% % lam is the linearized growth rate f'(0) of the field.
% %
% % Parameters lam, Dr, d, mu and nu can be set; default Dr = 2, d = 1,
% % mu = nu = 1, lam = 0.9*lam*.
% %
% % The paper gives a spreading-speed threshold (Dr = 2*d) [read], not a
% % stability threshold. Threshold [derived, finite differences]: sine
% % modes in s1, finite differences in s2 with a ghost node for the Robin
% % condition, modes n = 1..3, default Dr, d, mu, nu:
% %   lam* = 13.919655 (N = 40), 13.919824 (N = 60); critical mode n = 1.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), default Dr, d, mu, nu:
% %   - lam = 0.9*lam*: certified with the stock light settings (317,716
% %     decision variables) and with the proposed light row of the
% %     settings proposal of 11/15/2026 (319,545).
% %   - bisection on [0.5,1]*lam* (lam* at N = 60), 8 steps: both certify
% %     0.99805*lam* (lam = 13.892636), the top of the bracket.
% %   - lam = 1.02*lam*, 1.1*lam*: infeasible (stock light).
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

Dr = 2;     d = 1;      mu = 1;     nu = 1;
lam = 0.9*13.919824;
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
r = pde_var(s1,[0,1]);
f = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(r,'t')==Dr*diff(r,s1,2)+nu*subs(f,s2,0)-mu*r;
         diff(f,'t')==d*diff(f,s1,2)+d*diff(f,s2,2)+lam*f;
         subs(r,s1,0)==0;           subs(r,s1,1)==0;
         subs(f,s1,0)==0;           subs(f,s1,1)==0;
         subs(diff(f,s2),s2,0)==(nu/d)*subs(f,s2,0)-(mu/d)*r;
         subs(f,s2,1)==0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{berestycki2013influence,
%   title={The influence of a line with fast diffusion on Fisher-KPP propagation},
%   author={Berestycki, Henri and Roquejoffre, Jean-Michel and Rossi, Luca},
%   journal={Journal of Mathematical Biology},
%   volume={66},
%   pages={743--766},
%   year={2013},
%   note={arXiv:1210.3721; journal DOI 10.1007/s00285-012-0508-1 unverified}
% }
