function [PDE_t, st] = PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge(GUI,params) % MMP, 11/18/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=0.5;','sig=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) example-specific settings structure (header)     % MMP, 11/18/2026
%
% %---------------------------------------------------------------------% %
% % 2D heat equation with a dynamic boundary condition of reactive-
% % diffusive (Laplace-Beltrami) type on the edge s2 = 1 (box version of
% % Goldstein 2006; Vazquez, Vitillaro 2011) [derived]:
% % PDE         u_{t}(t,s1,s2) = u_{s1s1} + u_{s2s2} + lam*u;
% %                                               (s1,s2) in [0,1]^2
% %             v_{t}(t,s1)    = sig*v_{s1s1} - u_{s2}(t,s1,s2=1) + lam*v;
% %                                               s1 in [0,1]
% % With BCs    u_{s1} = 0 at s1 = 0, 1;   v_{s1} = 0 at s1 = 0, 1;
% %             u(t,s1,s2=0) = 0;          u(t,s1,s2=1) = v(t,s1).
% % The edge state v is the Dirichlet datum of u on s2 = 1, and the flux
% % of u drives v.
% %
% % Parameters lam and sig can be set; default sig = 1, lam = 0.9*lam*.
% %
% % Threshold [derived from the modes, exact]: in the cosine modes in s1,
% % mode n = 0 gives u = sin(k*s2), v = sin(k), with k*tan(k) = 1; modes
% % n >= 1 subtract (n*pi)^2 from every eigenvalue (sig = 1). Hence
% %   stable iff lam < lam* = k0^2 = 0.740174,   k0 = 0.860334.
% % Finite differences of the modes n = 0..3: 0.740247 (N = 40),
% % 0.740207 (N = 60). No threshold is published for a box.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), sig = 1:
% %   - lam = 0.9*lam*: certified with the stock light settings (426,044
% %     decision variables) and with the proposed light row of the
% %     settings proposal of 11/15/2026 (430,265).
% %   - bisection on [0.5,1]*lam*, 8 steps: both certify 0.99805*lam*
% %     (lam = 0.738728), the top of the bracket.
% %   - lam = 1.02*lam*, 1.1*lam*: UNKNOWN (stock light, after one raise).
% %                                                                         % MMP, 11/17/2026
% % No second output st: the light tier certifies this example        % MMP, 11/18/2026 (was); MMP, 11/17/2026
% % (record above). Its storage is the proposed light row since             % MMP, 11/17/2026
% % 11/17/2026; its slack of 11/17/2026 (wcross 0, 2-D faces at offset      % MMP, 11/17/2026
% % 1) was not rerun.                                                       % MMP, 11/17/2026
% %                                                                         % MMP, 11/18/2026
% % Example-specific settings of 11/18/2026 (output st; the paragraph       % MMP, 11/18/2026
% % above holds for the form 'P' only), for the direct form 'Q' that        % MMP, 11/18/2026
% % the library runs: lpisettings_sop('light') with the slack lift          % MMP, 11/18/2026
% % degree at its minimum, st.lpi.ineq.dD = 0. MOSEK, 11/18/2026, lam =     % MMP, 11/18/2026
% % 0.9*lam*: certified in 'Q' at 197,415 decision variables, 6,969         % MMP, 11/18/2026
% % rows, no raise (OPT_STATUS 0.996, eta 6.5e-10); light in 'Q' is         % MMP, 11/18/2026
% % refused by the size gate of the library pass (2,430,327 variables);     % MMP, 11/18/2026
% % light in 'P' certifies at 323,647 (OPT_STATUS 0.857); dD = 0 in 'P'     % MMP, 11/18/2026
% % is infeasible (295,135).                                                % MMP, 11/18/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/18/2026: Second output st for the direct form 'Q': light with
%                   the slack lift degree at its minimum (header);
% function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge(GUI,params) % MMP, 11/18/2026 (was)
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

sig = 1;    lam = 0.9*0.740174;
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
v = pde_var(s1,[0,1]);
u = pde_var([s1;s2],[0,1;0,1]);
PDE_t = [diff(v,'t')==sig*diff(v,s1,2)-subs(diff(u,s2),s2,1)+lam*v;
         diff(u,'t')==diff(u,s1,2)+diff(u,s2,2)+lam*u;
         subs(diff(v,s1),s1,0)==0;  subs(diff(v,s1),s1,1)==0;
         subs(diff(u,s1),s1,0)==0;  subs(diff(u,s1),s1,1)==0;
         subs(u,s2,0)==0;           subs(u,s2,1)==v];


% Example-specific settings (header): light, slack lift degree at its       % MMP, 11/18/2026
% minimum (dD = 0), for the form 'Q'                                        % MMP, 11/18/2026
if nargout>1                                                                % MMP, 11/18/2026
    st = lpisettings_sop('light');                                          % MMP, 11/18/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
end                                                                         % MMP, 11/18/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{goldstein2006derivation,
%   title={Derivation and physical interpretation of general boundary conditions},
%   author={Goldstein, Gis{\`e}le Ruiz},
%   journal={Advances in Differential Equations},
%   volume={11},
%   number={4},
%   pages={457--480},
%   year={2006}
% }
% @article{vazquez2011heat,
%   title={Heat equation with dynamical boundary conditions of reactive-diffusive type},
%   author={V{\'a}zquez, Juan Luis and Vitillaro, Enzo},
%   journal={Journal of Differential Equations},
%   volume={250},
%   pages={2143--2161},
%   year={2011},
%   doi={10.1016/j.jde.2010.12.012}
% }
% (Goldstein 2006: one repository lists pages 419-456; the DOI of
% Vazquez, Vitillaro 2011 is unverified.)
