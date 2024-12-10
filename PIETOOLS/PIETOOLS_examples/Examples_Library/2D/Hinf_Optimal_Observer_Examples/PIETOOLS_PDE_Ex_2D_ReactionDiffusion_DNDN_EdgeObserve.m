function [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
%
% %---------------------------------------------------------------------% %
% % Reaction-diffusion input-output-PDE with 
% %     Dirichlet-Neumann Dirichlet-Neumann BCs
% % PDE         x_{t}   = r*x + nu*x_{s1s1} + nu*x_{s2s2} + Cw*w1;
% % OUT         y1      = x(s2=d) + w2;
% %             y2      = x(s1=b) + w3;
% %             z       = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% % With BCs    x(s1=a) = 0;        x(s2=c) = 0;
% %             x_{s1}(s1=b) = 0;   x_{s2}(s2=d) = 0;
% %
% % Parameters nu and r, and limits a, b, c, d can be set.
% % Function Cw can also be set, defaults to
% %     Cw(s1,s2) = (s1^2-b)*(s2^2-d);
% % Stable for r <= mu, where
% %     mu = nu*pi^2*((1/4)/(b-a)^2 +(1/4)/(d-c)^2)
% %
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_estimator = 1;');

nu = 1;     r = 8;   ne = 1;
a = 0;  b=1;        c = 0;       d = 1;
Cw = (s1^2-b)*(s2^2-d);
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Term-based input format
PDE_t.x{1}.vars = [s1;s2];      PDE_t.x{1}.dom = [a,b;c,d];
PDE_t.x{1}.size = ne;
PDE_t.w{1}.vars = [];           PDE_t.w{1}.size = ne;
PDE_t.w{2}.vars = [];           PDE_t.w{2}.size = ne;
PDE_t.w{3}.vars = [];           PDE_t.w{3}.size = ne;
PDE_t.y{1}.vars = s1;           PDE_t.y{1}.dom = [a,b];
PDE_t.y{1}.size = ne;
PDE_t.y{2}.vars = s2;           PDE_t.y{2}.dom = [c,d];
PDE_t.y{2}.size = ne;
PDE_t.z{1}.vars = [];           PDE_t.z{1}.size = ne;

% PDE: x_{t} = [r, c1, c2] * [x; x_{s1s1}; x_{s2s2}] + (s1^2-1)*(s2^2-1)*w
PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
PDE_t.x{1}.term{1}.C = [r*eye(ne), nu*eye(ne), nu*eye(ne)];

PDE_t.x{1}.term{2}.w = 1;
PDE_t.x{1}.term{2}.C = Cw;

% OUT: y2(t,s2) = x(t,s1=b,s2) + w3;
PDE_t.y{2}.term{1}.x = 1;   PDE_t.y{2}.term{1}.loc = [b,s2];
PDE_t.y{2}.term{2}.w = 3;

% OUT: y1(t,s1) = x(t,s1,s2=d) + w2;
PDE_t.y{1}.term{1}.x = 1;   PDE_t.y{1}.term{1}.loc = [s1,d];
PDE_t.y{1}.term{2}.w = 2;

% OUT: z(t) = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
PDE_t.z{1}.term{1}.x = 1;
PDE_t.z{1}.term{1}.I = {[a,b];[c,d]};

% BC1: 0 = x(s1,c)                     % BC3: 0 = x_{s2}(s1,d)
PDE_t.BC{1}.term{1}.loc = [s1,c];      PDE_t.BC{3}.term{1}.loc = [s1,d];
                                       PDE_t.BC{3}.term{1}.D = [0,1];
% BC2: 0 = x(a,s2)                     % BC4: 0 = x_{s1}(b,s2)
PDE_t.BC{2}.term{1}.loc = [a,s2];      PDE_t.BC{4}.term{1}.loc = [b,s2];
                                       PDE_t.BC{4}.term{1}.D = [1,0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end