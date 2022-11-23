function [PDE_t] = PIETOOLS_PDE_Ex_2D_NS_Antonelli(GUI,params)
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
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
%
% %---------------------------------------------------------------------% %
% % Isentropic compressible Navier-Stokes, linearized around pE=1, and
% % vE=[vE1; vE2] = [s2;0] (Antonelli, 2021 [17])
% % PDE         p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}
% %             v1_{t} = -s2*v1_{s1} - v2 - (1/M^2)*p_{s1} + nu*(v1_{s1s1} + v1_{s2s2}) + lam*(v1_{s1s1} + v2_{s2s1})
% %             v2_{t} = -s2*v2_{s1}      - (1/M^2)*p_{s2} + nu*(v2_{s1s1} + v2_{s2s2}) + lam*(v1_{s1s2} + v2_{s2s2})
% % With BCs    p(s1=0) = 0;          p(s2=0) = 0;
% %             v(s1=0) = 0;          v(s2=0) = 0;
% %             v(s1=1) = 0;          v(s2=1) = 0;
% %
% % Parameters M and lam can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

M = 0.1;  lam = 1;  nu = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% Term-based input format
PDE_t.x{1}.vars = [s1;s2];   PDE_t.x{1}.dom = [0,1;0,1];
PDE_t.x{2}.vars = [s1;s2];   PDE_t.x{2}.dom = [0,1;0,1];
PDE_t.x{3}.vars = [s1;s2];   PDE_t.x{3}.dom = [0,1;0,1];

% PDE: x1_{t} = -s2*x1_{s1}
PDE_t.x{1}.term{1}.x = 1;
PDE_t.x{1}.term{1}.C = -s2;

% PDE: x1_{t} = ... - x2
PDE_t.x{1}.term{2}.x = 2;
PDE_t.x{1}.term{2}.C = -1;

% PDE: x1_{t} = ... -(1/M^2)*p_{s1}
PDE_t.x{1}.term{3}.x = 3;
PDE_t.x{1}.term{3}.C = -(1/M^2);

% PDE: x1_{t} = ... + nu*(v1_{s1s1} + v1_{s2s2})
PDE_t.x{1}.term{4}.x = [1; 1];
PDE_t.x{1}.term{4}.D = [2,0; 0,2];
PDE_t.x{1}.term{4}.C = [nu, nu];

% PDE: x1_{t} = ... + lam*(v1_{s1s2} + v2_{s2s2})
PDE_t.x{1}.term{5}.x = [1; 2];
PDE_t.x{1}.term{5}.D = [1,1; 0,2];
PDE_t.x{1}.term{5}.C = [lam, lam];


% PDE: x2_{t} = -s2*x2_{s1}
PDE_t.x{2}.term{1}.x = 2;
PDE_t.x{2}.term{1}.C = -s2;

% PDE: x2_{t} = ... -(1/M^2)*p_{s1}
PDE_t.x{2}.term{2}.x = 3;
PDE_t.x{2}.term{2}.C = -(1/M^2);

% PDE: x2_{t} = ... + nu*(v2_{s1s1} + v2_{s2s2})
PDE_t.x{2}.term{3}.x = [2; 2];
PDE_t.x{2}.term{3}.D = [2,0; 0,2];
PDE_t.x{2}.term{3}.C = [nu, nu];

% PDE: x2_{t} = ... + lam*(v1_{s1s1} + v2_{s2s1})
PDE_t.x{2}.term{4}.x = [1; 2];
PDE_t.x{2}.term{4}.D = [2,0; 1,1];
PDE_t.x{2}.term{4}.C = [lam, lam];


% PDE: x3_{t} = -s2*x3
PDE_t.x{3}.term{1}.x = 3;
PDE_t.x{3}.term{1}.C = -s2;

% PDE: x3_{t} = ... - x1_{s1} - x2_{s2}
PDE_t.x{3}.term{2}.x = [1;2];
PDE_t.x{3}.term{2}.D = [1,0; 0,1];
PDE_t.x{3}.term{2}.C = [-1, -1];


% BC1: 0 = x1(0,s2);                   % BC2: 0 = x1(1,s2);
PDE_t.BC{1}.term{1}.x = 1;             PDE_t.BC{2}.term{1}.x = 1;
PDE_t.BC{1}.term{1}.loc = [0,s2];      PDE_t.BC{2}.term{1}.loc = [1,s2];
% BC3: 0 = x1(s1,0);                   % BC4: 0 = x1(s1,1);
PDE_t.BC{3}.term{1}.x = 1;             PDE_t.BC{4}.term{1}.x = 1;
PDE_t.BC{3}.term{1}.loc = [s1,0];      PDE_t.BC{4}.term{1}.loc = [s1,1];

% BC5: 0 = x2(0,s2);                   % BC6: 0 = x2(1,s2);
PDE_t.BC{5}.term{1}.x = 2;             PDE_t.BC{6}.term{1}.x = 2;
PDE_t.BC{5}.term{1}.loc = [0,s2];      PDE_t.BC{6}.term{1}.loc = [1,s2];
% BC7: 0 = x2(s1,0);                   % BC8: 0 = x2(s1,1);
PDE_t.BC{7}.term{1}.x = 2;             PDE_t.BC{8}.term{1}.x = 2;
PDE_t.BC{7}.term{1}.loc = [s1,0];      PDE_t.BC{8}.term{1}.loc = [s1,1];

% BC5: 0 = x3(0,s2);                   % BC6: 0 = x3(s1,0);
PDE_t.BC{9}.term{1}.x = 3;             PDE_t.BC{10}.term{1}.x = 3;
PDE_t.BC{9}.term{1}.loc = [0,s2];      PDE_t.BC{10}.term{1}.loc = [s1,0];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{antonelli2021linear,
%  title={Linear stability analysis of the homogeneous Couette flow in a 2D isentropic compressible fluid},
%  author={Antonelli, Paolo and Dolce, Michele and Marcati, Pierangelo},
%  journal={Annals of PDE},
%  volume={7},
%  number={2},
%  pages={1--53},
%  year={2021},
%  publisher={Springer}
%}