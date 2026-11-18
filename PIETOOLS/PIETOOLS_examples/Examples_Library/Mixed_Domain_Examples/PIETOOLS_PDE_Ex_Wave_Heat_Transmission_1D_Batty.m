function [PDE_t] = PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings (this example has none).
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % Wave-heat transmission in 1D (Batty, Paunonen, Seifert 2016, Dirichlet
% % variant of Sec. 5; Zhang, Zuazua 2004) [read]:
% % PDE         u_{tt} = u_{xx}  on x in [-1,0];   w_{t} = w_{xx}  on x in [0,1]
% % With BCs    u(t,-1) = 0;   w(t,1) = 0;
% %             u_{t}(t,0) = w(t,0);   u_{x}(t,0) = w_{x}(t,0).
% %
% % Reflect the wave onto [0,1] (s1 = -x) and let e = u_{s1}, v = u_{t}:
% % PDE         e_{t} = v_{s1};   v_{t} = e_{s1};   w_{t} = w_{s1s1};
% %                                               s1 in [0,1]
% % With BCs    v(t,s1=1) = 0;   e(t,s1=0) = -w_{s1}(t,0);
% %             w(t,s1=1) = 0;   w(t,s1=0) = v(t,s1=0).
% %
% % Known result [read]: the energy is non-increasing and decays as
% % o(t^{-4}) for classical solutions, a sharp rate; the semigroup is not
% % exponentially stable (spectrum in the open left half-plane,
% % approaching the imaginary axis). No parameter threshold.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026):
% %   - epneg = 0 (non-increase of a weighted norm of Tx): certified with
% %     the stock light settings and with the proposed light row of the
% %     settings proposal of 11/15/2026 (2,520 decision variables each).
% %   - epneg = 0.01 and 0.001 (decay rate): infeasible (stock light,
% %     9,000 variables after two raises), as the polynomial decay
% %     requires.
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
pvar s1

%%% Executive Function:
evalin('base','stability = 1;');

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
e = pde_var(s1,[0,1],1);
v = pde_var(s1,[0,1],1);
w = pde_var(s1,[0,1]);
PDE_t = [diff(e,'t')==diff(v,s1);
         diff(v,'t')==diff(e,s1);
         diff(w,'t')==diff(w,s1,2);
         subs(v,s1,1)==0;           subs(e,s1,0)==-subs(diff(w,s1),s1,0);
         subs(w,s1,1)==0;           subs(w,s1,0)==subs(v,s1,0)];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{batty2016optimal,
%   title={Optimal energy decay in a one-dimensional coupled wave-heat system},
%   author={Batty, Charles J K and Paunonen, Lassi and Seifert, David},
%   journal={Journal of Evolution Equations},
%   volume={16},
%   pages={649--664},
%   year={2016},
%   note={arXiv:1505.04951; DOI 10.1007/s00028-015-0316-y unverified}
% }
% @article{zhang2004polynomial,
%   title={Polynomial decay and control of a 1-d hyperbolic-parabolic coupled system},
%   author={Zhang, Xu and Zuazua, Enrique},
%   journal={Journal of Differential Equations},
%   volume={204},
%   pages={380--438},
%   year={2004},
%   doi={10.1016/j.jde.2004.02.004}
% }
