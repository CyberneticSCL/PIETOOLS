close all; clc; clear all; 
% path(pathdef); 
% define pde parameters
E=68; I= 7.9e-5; rho=2.57; k=0.88; G=25; A =0.03141593;
alpha = E*I/(rho^2*I/k*G);
beta = rho*A/(rho^2*I/k*G);
gamma = (rho*I+E*I*rho/k*G)/(rho^2*I/k*G);

% define the PDE
pvar s theta; 

PDE_t.dom = [0,1];   PDE_t.n.n_pde = [1,0,2,0,1];

PDE_t.PDE.A{1}.coeff = [0 -beta]; PDE_t.PDE.A{1}.Lstate = 0; PDE_t.PDE.A{1}.Rstate = 2; PDE_t.PDE.A{1}.D = 0;
PDE_t.PDE.A{2}.coeff = [0;1]; PDE_t.PDE.A{2}.Lstate = 2; PDE_t.PDE.A{2}.Rstate = 0; PDE_t.PDE.A{2}.D = 0;
PDE_t.PDE.A{6}.coeff = [0 1;0 0]; PDE_t.PDE.A{6}.Lstate = 2; PDE_t.PDE.A{6}.Rstate = 2; PDE_t.PDE.A{6}.D = 0;
PDE_t.PDE.A{3}.coeff = [1 0]; PDE_t.PDE.A{3}.Lstate = 4; PDE_t.PDE.A{3}.Rstate = 2; PDE_t.PDE.A{3}.D = 0;
PDE_t.PDE.A{4}.coeff = [0 gamma]; PDE_t.PDE.A{4}.Lstate = 0; PDE_t.PDE.A{4}.Rstate = 2; PDE_t.PDE.A{4}.D = 2;
PDE_t.PDE.A{5}.coeff = -alpha; PDE_t.PDE.A{5}.Lstate = 0; PDE_t.PDE.A{5}.Rstate = 4; PDE_t.PDE.A{5}.D = 4;

% BCs
PDE_t.BC.Ebb{1}.coeff = [[1 0; 0 1]; zeros(6,2)]; PDE_t.BC.Ebb{1}.Rstate = 2; PDE_t.BC.Ebb{1}.delta = 0; PDE_t.BC.Ebb{1}.D = 0;
PDE_t.BC.Ebb{2}.coeff = [0; 0; 1; zeros(5,1)]; PDE_t.BC.Ebb{2}.Rstate = 4; PDE_t.BC.Ebb{2}.delta = 0; PDE_t.BC.Ebb{2}.D = 0;
PDE_t.BC.Ebb{3}.coeff = [zeros(3,2);eye(2);zeros(3,2)]; PDE_t.BC.Ebb{3}.Rstate = 2; PDE_t.BC.Ebb{3}.delta = 0;PDE_t.BC.Ebb{3}.D = 1;
PDE_t.BC.Ebb{4}.coeff = [zeros(5,1);1;0;0]; PDE_t.BC.Ebb{4}.Rstate = 4; PDE_t.BC.Ebb{4}.delta = 0;PDE_t.BC.Ebb{4}.D = 1;
PDE_t.BC.Ebb{5}.coeff = [zeros(6,1);1;0]; PDE_t.BC.Ebb{5}.Rstate = 4; PDE_t.BC.Ebb{5}.delta = 1;PDE_t.BC.Ebb{5}.D = 0;
PDE_t.BC.Ebb{6}.coeff = [zeros(7,1);1]; PDE_t.BC.Ebb{6}.Rstate = 4; PDE_t.BC.Ebb{6}.delta = 1;PDE_t.BC.Ebb{6}.D = 1;
PDE_t.BC.Ebb{7}.coeff = [zeros(6,1);-1;0]; PDE_t.BC.Ebb{7}.Rstate = 4; PDE_t.BC.Ebb{7}.delta = 1;PDE_t.BC.Ebb{7}.D = 2;
PDE_t.BC.Ebb{8}.coeff = [zeros(7,1);-1]; PDE_t.BC.Ebb{8}.Rstate = 4; PDE_t.BC.Ebb{8}.delta = 1;PDE_t.BC.Ebb{8}.D = 3;

PDE = PDE_t;

% convert PDE to PIE
use_pie=0;
if ~use_pie
if isfield(PDE,'n0')
    PIE = convert_PIETOOLS_PDE_batch(PDE);
elseif isfield(PDE,'n') && isfield(PDE.n,'n_pde')
    PIE = convert_PIETOOLS_PDE_terms(PDE);
else
    error('PDE not properly specified')
end
end

% rescale PIE system to [-1,1] interval
%PIE = rescalePIE(PIE,[-1,1]);

% define the simulation parameters
syms st sx;
n_pde = [1,0,2,0,1];
opts.N=8;
opts.dt=1e-3;
opts.tf=1;
opts.intScheme=1;
opts.Norder=2;
opts.Nsteps=opts.tf/opts.dt;
uinput.a = 0; uinput.b = 1; uinput.w = 0;
uinput.ic.ODE = 0; uinput.ic.PDE = [0;0;0;67*sx^2/13-42*sx^3/13+sx^4]; uinput.u = 0;
uinput.ifexact = 0;



% run simulation
sol = executive_PIESIM(PIE,opts,uinput,n_pde);