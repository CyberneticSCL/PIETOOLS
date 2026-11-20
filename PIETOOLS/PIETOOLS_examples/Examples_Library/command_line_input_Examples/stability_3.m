%   PDE: x_{t} = Mm*x - Lm*x_{s}                	| Mm, Lm, K00, K01,             Diagne 2012 [2] 
%   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]         | K10 and K11 may be
%        [x_-(s=1)] = [K10, K11] [x_-(s=0)]         | specified
%   The implementation is based on a Linearized Saint–Venant–Exner Model
%   The dissipativity condition of Diagne 2012 [2] fails for default
%   parameters  k1 = -.5;       k2 = -4.5;                                      
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable )
%
% The script declares the 3x3 hyperbolic system x_t = Mm*x - Lm*x_s on
% [0,1] with the boundary feedback [x1(1); x2(0); x3(0)] =
% K*[x1(0); x2(1); x3(1)], K = [K00, K01; K10, K11], where Lm, Mm and K are
% computed from the parameters Hs, Vs, Bs, g, Cf and the gains k1, k2;
% displays it, converts it to a PIE and displays the PIE. Workspace: PDE
% ('pde_struct') and PIE ('pie_struct'). If a cell array params of
% assignment strings (e.g. {'Vs=2;'}) exists in the workspace, its entries
% are evaluated after the default parameters Hs, Vs, Bs, g, Cf.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_3
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
% Initial coding SS - 11/16/2022
% MMP, 11/19/2026: Release comment pass (comment-only).

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
% Define dependent variables and system variable
x = pde_var(3,s,[0,1]); % x = [x_+; x_-];
%% Specify the parameters
 Hs = 1;       Vs = 1;     Bs = 0.1;
 g = 10;       Cf = 1;
 if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
 end
 % Characteristic speeds Lm (sorted), source matrix Mm, and the boundary
 % feedback K of the gains k1, k2
 Lm = diag(sort(eig([Vs,Hs,0;10,Vs,10;0,1*Vs^2,0])));
 lm1 = Lm(1,1);    lm2 = Lm(2,2);      lm3 = Lm(3,3);
 al1 = Cf*(3*Vs - 2*lm1) * Vs/Hs * lm1/((lm1-lm2)*(lm1-lm3));
 al2 = Cf*(3*Vs - 2*lm2) * Vs/Hs * lm2/((lm2-lm3)*(lm2-lm1));
 al3 = Cf*(3*Vs - 2*lm3) * Vs/Hs * lm3/((lm3-lm1)*(lm3-lm2));
 Mm = [al1*ones(3,1) , al2*ones(3,1) , al3*ones(3,1)];
 % Coefficients of the boundary feedback law
 a21 = (lm1 - lm2) * (1+ ((lm1-Vs) * (lm2-Vs))/(g*Hs));
 a32 = (lm2 - lm3) * (1+ ((lm2-Vs) * (lm3-Vs))/(g*Hs));
 a13 = (lm3 - lm1) * (1+ ((lm3-Vs) * (lm1-Vs))/(g*Hs));
 c21 = (lm3/g) * (lm1 - lm2);
 c32 = (lm1/g) * (lm2 - lm3);
 c13 = (lm2/g) * (lm3 - lm1);
 % Gains k1, k2, the boundary matrix K and its blocks
 k1 = -.5;       k2 = -4.5;
 pi2 = (a21 - c21*k1)/(a32 - c32*k1);
 pi3 = (a13 - c13*k1)/(a32 - c32*k1);
 chi2 = ((lm2 - Vs)/(lm1 - Vs)) * ((g + (lm2-Vs)*k2)/(g + (lm1-Vs)*k2));
 chi3 = ((lm3 - Vs)/(lm1 - Vs)) * ((g + (lm3-Vs)*k2)/(g + (lm1-Vs)*k2));
 K = [0,chi2,chi3; pi2,0,0; pi3,0,0];
 K00 = K(1,1);   K01 = K(1,2:3);
 K10 = K(2:3,1); K11 = K(2:3,2:3);
%% Define equations
eq_PDE = diff(x,t)==Mm*x-Lm*diff(x,s); % 	PDE: x_{t} = Mm*x - Lm*x_{s}
%   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)],
%            [x_-(s=1)]     [K10, K11] [x_-(s=0)]
% in the reference order lambda1 > lambda2 > 0 > lambda3
% Exemple order: lambda1 < 0 < lambda2 < lambda3
% Incoming traces:  [x1(1); x2(0); x3(0)]
% Outgoing traces:  [x1(0); x2(1); x3(1)]
B1 = [[1,0,0]*subs(x,s,1);
      [0,1,0;
       0,0,1]*subs(x,s,0)];
B2 = [[1,0,0]*subs(x,s,0);
      [0,1,0;
       0,0,1]*subs(x,s,1)];
eq_BC = B1 == [K00,K01; K10,K11]*B2;
%% addequations to pde system; set control inputs/observed inputs, if any
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: x_{t} = Mm*x - Lm*x_{s}                	| Mm, Lm, K00, K01,             Diagne 2012 [2] 
% %   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]         | K10 and K11 may be
% %        [x_-(s=1)] = [K10, K11] [x_-(s=0)]         | specified
% %   The implementation is based on a Linearized Saint–Venant–Exner Model
% x1 = pde_var(s,[0,1]);
% x2 = pde_var(s,[0,1]);
% x3 = pde_var(s,[0,1]);
% 
%  Hs = 1;       Vs = 1;     Bs = 0.1;
%  g = 10;       Cf = 1;
%  Lm = diag(sort(eig([Vs,Hs,0;10,Vs,10;0,1*Vs^2,0])));
%  lm1 = Lm(1,1);    lm2 = Lm(2,2);      lm3 = Lm(3,3);
%  al1 = Cf*(3*Vs - 2*lm1) * Vs/Hs * lm1/((lm1-lm2)*(lm1-lm3));
%  al2 = Cf*(3*Vs - 2*lm2) * Vs/Hs * lm2/((lm2-lm3)*(lm2-lm1));
%  al3 = Cf*(3*Vs - 2*lm3) * Vs/Hs * lm3/((lm3-lm1)*(lm3-lm2));
%  Mm = [al1*ones(3,1) , al2*ones(3,1) , al3*ones(3,1)];
%  a21 = (lm1 - lm2) * (1+ ((lm1-Vs) * (lm2-Vs))/(g*Hs));
%  a32 = (lm2 - lm3) * (1+ ((lm2-Vs) * (lm3-Vs))/(g*Hs));
%  a13 = (lm3 - lm1) * (1+ ((lm3-Vs) * (lm1-Vs))/(g*Hs));
%  c21 = (lm3/g) * (lm1 - lm2);
%  c32 = (lm1/g) * (lm2 - lm3);
%  c13 = (lm2/g) * (lm3 - lm1);
%  k1 = 1;       k2 = 2;
%  pi2 = (a21 - c21*k1)/(a32 - c32*k1);
%  pi3 = (a13 - c13*k1)/(a32 - c32*k1);
%  chi2 = ((lm2 - Vs)/(lm1 - Vs)) * ((g + (lm2-Vs)*k2)/(g + (lm1-Vs)*k2));
%  chi3 = ((lm3 - Vs)/(lm1 - Vs)) * ((g + (lm3-Vs)*k2)/(g + (lm1-Vs)*k2));
%  K = [0,chi2,chi3; pi2,0,0; pi3,0,0];
%  K00 = K(1,1);   K01 = K(1,2:3);
%  K10 = K(2:3,1); K11 = K(2:3,2:3);
% 
% %% Define equations
% eq_PDE = [diff(x1,t)==Mm(1,:)*[x1;x2;x3]-Lm(1,:)*diff([x1;x2;x3],s);
%           diff(x2,t)==Mm(2,:)*[x1;x2;x3]-Lm(2,:)*diff([x1;x2;x3],s);
%           diff(x3,t)==Mm(3,:)*[x1;x2;x3]-Lm(3,:)*diff([x1;x2;x3],s)]; % 	PDE: x_{t} = Mm*x - Lm*x_{s}
% eq_BC = [subs(x1,s,0); subs([x2;x3],s,1)]==[K00,K01;K10,K11]*[subs([x1;x2],s,1);subs(x3,s,0)];    
% %   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)], [x_-(s=1)] = [K10, K11] [x_-(s=0)]
% %% initialize pde system;
% pde = initialize([eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);



% %% Define dependent variables and system variable
% %   PDE: x_{t} = Mm*x - Lm*x_{s}                	| Mm, Lm, K00, K01,             Diagne 2012 [2] 
% %   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]         | K10 and K11 may be
% %        [x_-(s=1)] = [K10, K11] [x_-(s=0)]         | specified
% %   The implementation is based on a Linearized Saint–Venant–Exner Model
% x = state('pde',3); % x = [x_+; x_-];
% pde = sys();
% 
%  Hs = 1;       Vs = 1;     Bs = 0.1;
%  g = 10;       Cf = 1;
%  Lm = diag(sort(eig([Vs,Hs,0;10,Vs,10;0,1*Vs^2,0])));
%  lm1 = Lm(1,1);    lm2 = Lm(2,2);      lm3 = Lm(3,3);
%  al1 = Cf*(3*Vs - 2*lm1) * Vs/Hs * lm1/((lm1-lm2)*(lm1-lm3));
%  al2 = Cf*(3*Vs - 2*lm2) * Vs/Hs * lm2/((lm2-lm3)*(lm2-lm1));
%  al3 = Cf*(3*Vs - 2*lm3) * Vs/Hs * lm3/((lm3-lm1)*(lm3-lm2));
%  Mm = [al1*ones(3,1) , al2*ones(3,1) , al3*ones(3,1)];
%  a21 = (lm1 - lm2) * (1+ ((lm1-Vs) * (lm2-Vs))/(g*Hs));
%  a32 = (lm2 - lm3) * (1+ ((lm2-Vs) * (lm3-Vs))/(g*Hs));
%  a13 = (lm3 - lm1) * (1+ ((lm3-Vs) * (lm1-Vs))/(g*Hs));
%  c21 = (lm3/g) * (lm1 - lm2);
%  c32 = (lm1/g) * (lm2 - lm3);
%  c13 = (lm2/g) * (lm3 - lm1);
%  k1 = 1;       k2 = 2;
%  pi2 = (a21 - c21*k1)/(a32 - c32*k1);
%  pi3 = (a13 - c13*k1)/(a32 - c32*k1);
%  chi2 = ((lm2 - Vs)/(lm1 - Vs)) * ((g + (lm2-Vs)*k2)/(g + (lm1-Vs)*k2));
%  chi3 = ((lm3 - Vs)/(lm1 - Vs)) * ((g + (lm3-Vs)*k2)/(g + (lm1-Vs)*k2));
%  K = [0,chi2,chi3; pi2,0,0; pi3,0,0];
%  K00 = K(1,1);   K01 = K(1,2:3);
%  K10 = K(2:3,1); K11 = K(2:3,2:3);
% 
% %% Define equations
% eq_PDE = diff(x,t)==Mm*x-Lm*diff(x,s); % 	PDE: x_{t} = Mm*x - Lm*x_{s}
% eq_BC = [[1,0,0]*subs(x,s,0); [0,1,0;0,0,1]*subs(x,s,1)]==[K00,K01;K10,K11]*[[1,0,0;0,1,0]*subs(x,s,1);[0,0,1]*subs(x,s,0)];    
% %   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)], [x_-(s=1)] = [K10, K11] [x_-(s=0)]
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
