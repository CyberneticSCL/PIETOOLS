function [PDE_t] = PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN(GUI,params)
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
% % PDE         x_{t}   = r*x + nu*x_{s1s1} + nu*x_{s2s2} + w;
% % OUT         z       = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% % With BCs    x(s1=a) = 0;        x(s2=c) = 0;
% %             x_{s1}(s1=b) = 0;   x_{s2}(s2=d) = 0;
% %
% % Parameters nu and r, and limits a, b, c, d can be set.
% % Stable for r <= mu_{1,1}, where
% %     mu_{m,n} = nu*pi^2*((m-1/2)^2/(b-a)^2 +(n-1/2)^2/(d-c)^2)
% % Then, the L2-gain given by
% % gam = 8*sqrt((b-a)*(d-c))/pi^2 
% %             *sqrt(sum_{m=1}^{infty}sum_{n=1}^{infty} 1/((mu_{2m-1,2n-1}-r)*(2m-1)*(2n-1))^2
% %
% % For nu = 1, r = 1, [a,b]=[c,d]=[0,1], we have gam = 0.20668;
% %
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_gain = 1;');
%evalin('base','Hinf_gain_dual = 1;');

nu = 1;     r = 1;   ne = 1;
a = 0;  b = 1;       c = 0;  d = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input format
clear stateNameGenerator
x = pde_var(ne,[s1;s2],[a,b;c,d]);
w = pde_var('in',ne,[s1;s2],[a,b;c,d]);
z = pde_var('out',ne);
PDE_t = [diff(x,'t')==r*x+nu*(diff(x,s1,2)+diff(x,s2,2))+w;
         z==int(x,[s1;s2],[a,b;c,d]);
         subs(x,s1,a)==0;       subs(diff(x,s1),s1,b)==0;
         subs(x,s2,c)==0;       subs(diff(x,s2),s2,d)==0];


% %%% Term-based input format
% PDE_t.x{1}.vars = [s1;s2];      PDE_t.x{1}.dom = [a,b;c,d];
% PDE_t.x{1}.size = ne;
% PDE_t.w{1}.vars = [s1;s2];      PDE_t.w{1}.dom = [a,b;c,d];
% PDE_t.w{1}.size = ne;
% PDE_t.z{1}.vars = [];           PDE_t.z{1}.size = ne;
% 
% % PDE: x_{t} = [r, c1, c2] * [x; x_{s1s1}; x_{s2s2}] + w
% PDE_t.x{1}.term{1}.D = [0,0; 2,0; 0,2];
% PDE_t.x{1}.term{1}.C = [r*eye(ne), nu*eye(ne), nu*eye(ne)];
% 
% PDE_t.x{1}.term{2}.w = 1;
% 
% % OUT: z(t) = int_{a}^{b} int_{c}^{d} x(t,s1,s2) ds2 ds1
% PDE_t.z{1}.term{1}.x = 1;
% PDE_t.z{1}.term{1}.I = {[a,b];[c,d]};
% 
% % BC1: 0 = x(s1,c)                     % BC3: 0 = x_{s2}(s1,d)
% PDE_t.BC{1}.term{1}.loc = [s1,c];      PDE_t.BC{3}.term{1}.loc = [s1,d];
%                                        PDE_t.BC{3}.term{1}.D = [0,1];
% % BC2: 0 = x(a,s2)                     % BC4: 0 = x_{s1}(b,s2)
% PDE_t.BC{2}.term{1}.loc = [a,s2];      PDE_t.BC{4}.term{1}.loc = [b,s2];
%                                        PDE_t.BC{4}.term{1}.D = [1,0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

%%
% For comparison, compute the L2-gain (approximately)
M = 10;     N = 10;
Lx = b-a;   Ly = d-c;
r_lim = nu*pi^2*((1/4)/Lx^2 +(1/4)/Ly^2);       % Stable iff r<=r_lim
mu_mn = nu*pi^2*((((1:M)'-(1/2))/Lx).^2 +(((1:N)-(1/2))/Ly).^2);
mn_fact = ((2*(1:M)'-1)).*((2*(1:N)-1));
full_fact = 1./((mu_mn-r).*mn_fact);
fct = sqrt(sum(full_fact.^2,'all'));
L2_gain = (8*sqrt(Lx*Ly)/(pi^2))*fct;       % L2 gain (approximately)

end