function [PDE_t,PDE_b,st] = PIETOOLS_PDE_Ex_Diagne(GUI,params)              % DB, MMP, 11/21/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,PDE_b,ST] = PIETOOLS_PDE_EX_DIAGNE(GUI,PARAMS) returns the 3-state
% linear hyperbolic system of balance laws below (linearized
% Saint-Venant-Exner model with boundary feedback), in the term-based and
% in the batch format, and the recommended settings ST of its stability
% test. It sets stability = 1 and stability_form = 'P' (storage
% <T x, P T x>, P >= eppos I) in the base workspace. With GUI = 1 it also
% loads the associated GUI save file into a new PIETOOLS_PDE_GUI window.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % DB, MMP, 11/21/2026
%               for the tests this file names                               % DB, MMP, 11/21/2026
%
% %---------------------------------------------------------------------% %
% % Diagne 2012 example (see reference below): 
% % PDE          x_{t} = Mm*x - Lm*x_{s}
% % with BC [x_+(s=0)] = [K00, K01] [x_+(s=1)]
% %         [x_-(s=1)] = [K10, K11] [x_-(s=0)]
% %
% % The implementation is based on a Linearized Saint–Venant–Exner Model
% % Mm, Lm, K00, K01, K10 and K11 may be specified.
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_Diagne
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding DJ - 08/18/2022
% DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the BCs read
%   [x1(1); x2(0); x3(0)] = K [x1(0); x2(1); x3(1)], the incoming traces
%   in terms of the outgoing ones, in both formats (were [x1(0); x2(1);
%   x3(1)] = K [x1(1); x2(0); x3(0)]); Lm uses g in place of the constant
%   10; default gains k1 = -0.5, k2 = -4.5 (were 1, 2); new output st, the
%   recommended settings of the stability test; stability_form = 'P' set
%   in the base workspace; help and comments revised.
% function [PDE_t,PDE_b] = PIETOOLS_PDE_Ex_Diagne(GUI,params)               % DB, MMP, 11/21/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','stability = 1;')
evalin('base','stability_form = ''P'';');  % form of the stability test     % DB, MMP, 11/21/2026
%evalin('base','stability_dual = 1;')

npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Default Parameters
Hs = 1;       Vs = 1;     Bs = 0.1;
g = 10;       Cf = 1;

 % Characteristic speeds Lm and source matrix Mm of the linearized model
 Lm = diag(sort(eig([Vs,Hs,0;g,Vs,g;0,1*Vs^2,0])));
 lm1 = Lm(1,1);    lm2 = Lm(2,2);      lm3 = Lm(3,3);

 al1 = Cf*(3*Vs - 2*lm1) * Vs/Hs * lm1/((lm1-lm2)*(lm1-lm3));
 al2 = Cf*(3*Vs - 2*lm2) * Vs/Hs * lm2/((lm2-lm3)*(lm2-lm1));
 al3 = Cf*(3*Vs - 2*lm3) * Vs/Hs * lm3/((lm3-lm1)*(lm3-lm2));

 Mm = [al1*ones(3,1) , al2*ones(3,1) , al3*ones(3,1)];

 a21 = (lm1 - lm2) * (1+ ((lm1-Vs) * (lm2-Vs))/(g*Hs));
 a32 = (lm2 - lm3) * (1+ ((lm2-Vs) * (lm3-Vs))/(g*Hs));
 a13 = (lm3 - lm1) * (1+ ((lm3-Vs) * (lm1-Vs))/(g*Hs));
 c21 = (lm3/g) * (lm1 - lm2);
 c32 = (lm1/g) * (lm2 - lm3);
 c13 = (lm2/g) * (lm3 - lm1);

 % Boundary feedback with gains k1, k2
 % k1 = 1;       k2 = 2;
 k1 = -.5;       k2 = -4.5;
 pi2 = (a21 - c21*k1)/(a32 - c32*k1);
 pi3 = (a13 - c13*k1)/(a32 - c32*k1);
 chi2 = ((lm2 - Vs)/(lm1 - Vs)) * ((g + (lm2-Vs)*k2)/(g + (lm1-Vs)*k2));
 chi3 = ((lm3 - Vs)/(lm1 - Vs)) * ((g + (lm3-Vs)*k2)/(g + (lm1-Vs)*k2));

 K = [0,chi2,chi3; pi2,0,0; pi3,0,0];
 % K is built in the example ordering
 % [x_-(1); x_+(0)] = K [x_-(0); x_+(1)]
K11 = K(1,1);
K10 = K(1,2:3);
K01 = K(2:3,1);
K00 = K(2:3,2:3);
 % K00 = K(1,1);   K01 = K(1,2:3);
 % K10 = K(2:3,1); K11 = K(2:3,2:3);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
 
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

ne = size(Lm,1);


%%% Batch input format
PDE_b.n0 = 0;   PDE_b.n1 = ne;   PDE_b.n2 = 0;
PDE_b.dom = [0,1];

PDE_b.A0 = Mm;
PDE_b.A1 = -Lm;
% [x_-(1); x_+(0)] = K [x_-(0); x_+(1)]

ny1 = size(K00,1);   % number of x_+ states = 2
ny2 = size(K11,1);   % number of x_- states = 1

on0 = eye(ny1);
on1 = eye(ny2);
zer12 = zeros(ny1,ny2);

% Converts raw state [x_-; x_+] to paper order [x_+; x_-]
Ppm = [0 1 0;
       0 0 1;
       1 0 0];
C0_pm = [ on0,    -K01;
               zer12', -K11 ];

C1_pm = [ -K00,  zer12;
             -K10,  on1 ];

C0_raw = C0_pm * Ppm;
C1_raw = C1_pm * Ppm;
%ny1 = size(K00,1);   ny2 = size(K10,1);
%on0 = eye(ny1);   on1 = eye(ny2);   zer12 = zeros(ny1,ny2);
PDE_b.B = [C0_raw, C1_raw];


%%% Term-based input format
PDE_t.x{1}.vars = s;
PDE_t.x{1}.dom = [0,1];

% PDE: x_{t} = Mm * x
PDE_t.x{1}.term{1}.C = Mm;

% PDE: x_{t} = ... -Lm * x_{s}
PDE_t.x{1}.term{2}.D = 1;
PDE_t.x{1}.term{2}.C = -Lm;

% BCs: 0 = [1,-K01;0,-K11] * x(0)
PDE_t.BC{1}.term{1}.loc = 0;
PDE_t.BC{1}.term{1}.C = C0_raw;

% BCs: 0 = ... + [-K00,0;-K10,1] * x(1)
PDE_t.BC{1}.term{2}.loc = 1;
PDE_t.BC{1}.term{2}.C = C1_raw;

%%% Command-line format
% PDE_c=sys();
% x = state('pde',3);
% eq_PDE = diff(x,t)==Mm*x -Lm*diff(x,s); 
% eq_BC =  [on0, -K01; zer12', -K11]*subs(x,s,0) + [-K00, zer12; -K10, on1]*subs(x,s,1)==0;        
% % addequations to pde system; set control inputs/observed inputs, if any
% PDE_c = addequation(PDE_c,[eq_PDE;eq_BC]);
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(fullfile(root,'PIETOOLS_PDE_Ex_Diagne_GUI.mat'));
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

% Example-specific settings                                                 % DB, MMP, 11/21/2026
% (stripped tier; 1-D storage without the upper-integral basis operator,
% s-length 0, and theta-lengths 1 (multiplier), 1 (lower integral); slack
% lift degree at its minimum Dmin)
st = lpisettings_sop('stripped');                                           % DB, MMP, 11/21/2026
st.lpi.pos.len_com{1} = 0;                                                  % DB, MMP, 11/21/2026
st.lpi.pos.len_sos{1} = [1;2;1];                                            % DB, MMP, 11/21/2026
st.lpi.ineq.dD = 0;                                                         % DB, MMP, 11/21/2026

end
% @article{diagne2012lyapunov,
%   title={Lyapunov exponential stability of 1-D linear hyperbolic systems of balance laws},
%   author={Diagne, Ababacar and Bastin, Georges and Coron, Jean-Michel},
%   journal={Automatica},
%   volume={48},
%   number={1},
%   pages={109--114},
%   year={2012},
%   publisher={Elsevier}
% }