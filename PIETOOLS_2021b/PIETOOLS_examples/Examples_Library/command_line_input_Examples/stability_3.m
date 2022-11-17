clc; clear;
pvar s t theta; % define independent variables
%% Define dependent variables and system variable
%   PDE: x_{t} = Mm*x - Lm*x_{s}                	| Mm, Lm, K00, K01,             Diagne 2012 [2] 
%   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)]         | K10 and K11 may be
%        [x_-(s=1)] = [K10, K11] [x_-(s=0)]         | specified
%   The implementation is based on a Linearized Saint–Venant–Exner Model
x = state('pde',3); % x = [x_+; x_-];
pde = sys();

 Hs = 1;       Vs = 1;     Bs = 0.1;
 g = 10;       Cf = 1;
 Lm = diag(sort(eig([Vs,Hs,0;10,Vs,10;0,1*Vs^2,0])));
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
 k1 = 1;       k2 = 2;
 pi2 = (a21 - c21*k1)/(a32 - c32*k1);
 pi3 = (a13 - c13*k1)/(a32 - c32*k1);
 chi2 = ((lm2 - Vs)/(lm1 - Vs)) * ((g + (lm2-Vs)*k2)/(g + (lm1-Vs)*k2));
 chi3 = ((lm3 - Vs)/(lm1 - Vs)) * ((g + (lm3-Vs)*k2)/(g + (lm1-Vs)*k2));
 K = [0,chi2,chi3; pi2,0,0; pi3,0,0];
 K00 = K(1,1);   K01 = K(1,2:3);
 K10 = K(2:3,1); K11 = K(2:3,2:3);

%% Define equations
eq_PDE = diff(x,t)==Mm*x-Lm*diff(x,s); % 	PDE: x_{t} = Mm*x - Lm*x_{s}
eq_BC = [[1,0,0]*subs(x,s,0); [0,1,0;0,0,1]*subs(x,s,1)]==[K00,K01;K10,K11]*[[1,0,0;0,1,0]*subs(x,s,1);[0,0,1]*subs(x,s,0)];    
%   BCs: [x_+(s=0)] = [K00, K01] [x_+(s=1)], [x_-(s=1)] = [K10, K11] [x_-(s=0)]
%% addequations to pde system; set control inputs/observed inputs, if any
pde = addequation(pde,[eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);
