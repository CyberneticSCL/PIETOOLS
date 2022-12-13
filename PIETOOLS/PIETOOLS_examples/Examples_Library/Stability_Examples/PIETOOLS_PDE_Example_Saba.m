function [PDE_t,PDE_b] = PIETOOLS_PDE_Example_Saba(index,GUI,params)
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
% % Saba 2019 examples (see reference below): 
% % PDE        x1_{t} = sig1*x2 - (1/r1)*x1_{s} 
% %            x2_{t} = sig2*x1 + (1/r2)*x2_{s}
% % with BCs   x1(s=0) = qb*x2(s=0)
% %            x2(s=1) = pb*x1(s=1)
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s theta s1 s2 theta1 theta2

%%% Executive Function:
evalin('base','stability = 1;');
% evalin('base','stability_dual = 1;');


%%%% There are several examples included here. Add a decimal to
%%%% your example input to specify a particular set of parameters.
if index==4 || index==4.1
 r1=.8;   r2=1.1;   sig1=2.3;   sig2=-3.5;    qb=-.7;   pb=.5;  % Stable with stripped settings
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_1.mat');
 
elseif index==4.2
 r1=.5;   r2=1.1;   sig1=1;     sig2=-.1;     qb=1.2;   pb=0;   % Stable with stripped settings
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_2.mat');
 
elseif index==4.3
 r1=.5;   r2=1.1;   sig1=1;     sig2=1;       qb=1.2;   pb=-.7;
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_3.mat');
 
elseif index==4.4
 r1=.5;   r2=1.1;   sig1=1;     sig2=.663;    qb=1.2;   pb=0;   % max sig2=.663 
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_4.mat');
 
elseif index==4.5
 r1=.5;   r2=1.1;   sig1=1;     sig2=2.048;   qb=1.2;   pb=-.4; % max sig2=1.049 
 dir = fullfile(root,'Examples_PDE_GUI/Stability_Examples/Saba_Ex_5.mat');
 
else
 error(['Example ',num2str(indx),'does not exist, please choose from 4.1, 4.2, 4.3, 4.4 or 4.5']);
end

npars = length(params);
if npars~=0
 %%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end


 %%% Batch input format
 PDE_b.n0 = 0;   PDE_b.n1 = 2;   PDE_b.n2 = 0;
 PDE_b.dom = [0,1];

 PDE_b.A0 = [0 sig1; sig2 0]; 
 PDE_b.A1 = [-1/r1 0; 0 1/r2];
 PDE_b.B = [1 -qb, 0, 0; 0, 0, -pb, 1];

 
 %%% Term-based input format
 PDE_t.x{1}.vars = s;
 PDE_t.x{1}.dom = [0,1];
 
 % PDE: x_{t} = [0,sig1;sig2,0] * x
 PDE_t.x{1}.term{1}.C = [0, sig1; sig2, 0];

 % PDE: x_{t} = ... + [-1/r1,0;0,1/r2] * x_{s}
 PDE_t.x{1}.term{2}.D = 1;
 PDE_t.x{1}.term{2}.C = [-1/r1, 0; 0, 1/r2];
 
 % BCs: 0 = [1,-qb;0,0] * x(0)
 PDE_t.BC{1}.term{1}.loc = 0;
 PDE_t.BC{1}.term{1}.C = [1 -qb; 0 0];

 % BCs: 0 = ... + [0,0;-pb,1] * x(1)
 PDE_t.BC{1}.term{2}.loc = 1;
 PDE_t.BC{1}.term{2}.C = [0 0; -pb 1];

 
if GUI
    %%% Associated GUI save file
    app = PIETOOLS_PDE_GUI;
    load(dir);
    logval = app.loadData(data);
    if logval
        disp("Failed to load data object. Incorrect structure");
    end
end

end
% @article{saba2019stability,
%   title={Stability Analysis for a Class of Linear 2x2 Hyperbolic PDEs Using a Backstepping Transform},
%   author={Saba, David Bou and Argomedo, Federico Bribiesca and Auriol, Jean and Di Loreto, Michael and Di Meglio, Florent},
%   journal={IEEE Transactions on Automatic Control},
%   year={2019},
%   publisher={IEEE}
% }