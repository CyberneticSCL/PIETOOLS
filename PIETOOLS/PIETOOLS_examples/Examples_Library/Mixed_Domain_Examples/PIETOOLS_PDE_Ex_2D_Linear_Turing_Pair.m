function [PDE_t] = PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'d=3;','bc=''D'';'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % Linear Turing system: two 2D states coupled in-domain (Turing 1952;
% % Murray 2003, ch. 2), x = (u,v):
% % PDE         x_{t}(t,s1,s2) = diag(1,d)*(x_{s1s1} + x_{s2s2}) + J*x;
% %                     J = [20, -25; 30, -30];    (s1,s2) in [0,1]^2
% % With BCs    bc = 'N': x_{s1} = 0 at s1 = 0, 1;  x_{s2} = 0 at s2 = 0, 1;
% %             bc = 'D': x = 0 on all edges, with J replaced by 2*J.
% %
% % Parameters d and bc can be set; default bc = 'N', d = 0.9*d*.
% %
% % Threshold [derived from the modes, exact]: cosine mode with
% % k^2 = pi^2*(m^2+n^2) is stable iff tr_k = -10-(1+d)*k^2 < 0 and
% % det_k = (20-k^2)*(-30-d*k^2) + 750 > 0. tr_k < 0 for every k, and
% % det_k changes sign first at k^2 = pi^2:
% %   stable iff d < d* = 4.461640   (same d* for bc = 'D').
% % Maximum real part over the modes m,n <= 20: -0.759 (0.9*d*),
% % -0.070 (0.99*d*), +0.069 (1.01*d*).
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), bc = 'N', d = 0.9*d*: certified with the stock
% % light settings (2,102,410 decision variables, 740 s); the run with
% % the proposed light row of the settings proposal of 11/15/2026
% % (2,120,516 variables) was not completed.
% %                                                                         % MMP, 11/17/2026
% % No second output st: the light tier certified this example with the     % MMP, 11/17/2026
% % stock storage (record above); the light tier of 11/17/2026 (the         % MMP, 11/17/2026
% % proposed row, wcross 0, 2-D faces at offset 1) was not rerun.           % MMP, 11/17/2026
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

d = 0.9*4.461640;   bc = 'N';
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
J = [20 -25;30 -30];    Dm = [1 0;0 d];
x = pde_var(2,[s1;s2],[0,1;0,1]);
if bc=='N'
    PDE_t = [diff(x,'t')==Dm*diff(x,s1,2)+Dm*diff(x,s2,2)+J*x;
             subs(diff(x,s1),s1,0)==[0;0];  subs(diff(x,s1),s1,1)==[0;0];
             subs(diff(x,s2),s2,0)==[0;0];  subs(diff(x,s2),s2,1)==[0;0]];
else
    PDE_t = [diff(x,'t')==Dm*diff(x,s1,2)+Dm*diff(x,s2,2)+2*J*x;
             subs(x,s1,0)==[0;0];   subs(x,s1,1)==[0;0];
             subs(x,s2,0)==[0;0];   subs(x,s2,1)==[0;0]];
end


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{turing1952chemical,
%   title={The chemical basis of morphogenesis},
%   author={Turing, Alan Mathison},
%   journal={Philosophical Transactions of the Royal Society of London. Series B},
%   volume={237},
%   pages={37--72},
%   year={1952},
%   doi={10.1098/rstb.1952.0012}
% }
% @book{murray2003mathematical,
%   title={Mathematical Biology II: Spatial Models and Biomedical Applications},
%   author={Murray, James D},
%   edition={3},
%   publisher={Springer},
%   year={2003}
% }
% (DOI of Turing 1952 and the chapter of Murray 2003 unverified.)
