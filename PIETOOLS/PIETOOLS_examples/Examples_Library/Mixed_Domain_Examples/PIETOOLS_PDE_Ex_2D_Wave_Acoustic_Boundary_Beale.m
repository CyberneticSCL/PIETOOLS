function [PDE_t] = PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale(GUI,params)
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
% % 2D wave equation with an acoustic boundary condition on the edge
% % s2 = 0 (Beale 1976; Munoz Rivera, Qin 2003) [unverified standard form],
% % m = d = k = rho = c = 1:
% % PDE         phi_{tt} = phi_{s1s1} + phi_{s2s2};    (s1,s2) in [0,1]^2
% %             del_{tt} + del_{t} + del = -phi_{t}(t,s1,0);   s1 in [0,1]
% % With BCs    del_{t} = d phi/d nu = -phi_{s2}(t,s1,0);
% %             phi = 0 at s1 = 0, 1 and at s2 = 1.
% %
% % Let psi = phi_{t}, eta = del_{t}:
% % PDE         phi_{t} = psi;   psi_{t} = phi_{s1s1} + phi_{s2s2};
% %             del_{t} = eta;   eta_{t} = -eta - del - psi(t,s1,s2=0);
% % With BCs    phi = psi = 0 at s1 = 0, 1;
% %             phi_{s2}(t,s1,0) = -eta(t,s1);
% %             phi(t,s1,1) = 0;   psi(t,s1,1) = 0;
% %             del = eta = 0 at s1 = 0, 1.
% % Modification: the edge states del, eta are 2nd order in s1, with
% % del = eta = 0 at s1 = 0, 1 (compatible with phi = 0 there). The model
% % has no s1-derivative on the edge, but the converter rejects edge
% % states of order 0, since eta is the Neumann datum of phi.
% %
% % Known result [abstract]: the generator has essential spectrum and the
% % system is not uniformly (exponentially) stable (Beale 1976); the
% % energy decays polynomially on smooth domains (Munoz Rivera, Qin
% % 2003). No parameter threshold.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), epneg = 0: infeasible with the stock light
% % settings (2,032,912 decision variables after one raise, 737 s); the
% % run with the proposed light row of the settings proposal of
% % 11/15/2026 was not completed.
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

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
ph = pde_var([s1;s2],[0,1;0,1],[2;2]);
ps = pde_var([s1;s2],[0,1;0,1],[2;1]);          % 2nd order in s1, 1st in s2
de = pde_var(s1,[0,1],2);                       % order 2: edge datum of ph must be C2 in s1
et = pde_var(s1,[0,1],2);
PDE_t = [diff(ph,'t')==ps;
         diff(ps,'t')==diff(ph,s1,2)+diff(ph,s2,2);
         diff(de,'t')==et;
         diff(et,'t')==-et-de-subs(ps,s2,0);
         subs([ph;ps],s1,0)==0;     subs([ph;ps],s1,1)==0;
         subs(diff(ph,s2),s2,0)==-et;
         subs(ph,s2,1)==0;          subs(ps,s2,1)==0;
         subs([de;et],s1,0)==0;     subs([de;et],s1,1)==0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{beale1976spectral,
%   title={Spectral properties of an acoustic boundary condition},
%   author={Beale, J Thomas},
%   journal={Indiana University Mathematics Journal},
%   volume={25},
%   number={9},
%   pages={895--917},
%   year={1976},
%   doi={10.1512/iumj.1976.25.25071}
% }
% @article{munozrivera2003polynomial,
%   title={Polynomial decay for the energy with an acoustic boundary condition},
%   author={Mu{\~n}oz Rivera, Jaime E and Qin, Yuming},
%   journal={Applied Mathematics Letters},
%   volume={16},
%   number={2},
%   pages={249--256},
%   year={2003},
%   doi={10.1016/S0893-9659(03)80039-3}
% }
% (DOI of Beale 1976 unverified.)
