%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_PDE_library_PIETOOLS.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = examples_PDE_library_PIETOOLS_2D(varargin)
% This library file contains the definition of some common 2 DODE-PDE systems
% To use one of the examples, call this function with as FIRST input the 
% index of the example you wish to use. Please check the table of contents 
% below to see which index corresponds to which function. The output of the 
% function then corresponds to a Matlab structure describing the desired PDE.
%
% Inputs:
%     (output1,output2) = examples_PDE_library_PIETOOLS(example,option1,option2,option3,...)
%       - example - integer (or possibly double) number (1-27) corresponding to selected example problem    
%       - option1,option2,option3, etc.
%            - option# - 'batch' or 'BATCH' - output# will be the PDE data structure in batch format  
%            - option# - 'terms' or 'TERMS' - output# will be the PDE data structure in terms format
%            - option# - 'gui' or 'GUI' - no output. Executes PIETOOLS PDE GUI
%             and loads saved problem data into GUI.
%            - option# - 'lambda=XXX' overrides the default parameter value
%             lambda with value XXX for the given problem number. The list of allowable parameter 
%             overrides for each example problem is given in the table of
%             contents below. DO NOT ABUSE THIS POWER!
%             NOTE: Changing parameters does not change the GUI file!
%            
%
%
% Several examples come with parameters that may be adjusted. Check the
% "Parameters" column in the table of contents below. If you wish to adjust
% one of the parameter values, specify this adjustment using a character 
% array as one of your arguments. For example, add argument "lam = 10" to
% set the parameter "lam" equal to 10. 
% Also, most examples in each problem type include the option to solve the 
% LPI according to their problem type. If you want to invoke a different
% LPI, specify this option after running the function (e.g. set
% "stability = 1", "Hinf_gain=1", etc.).
%
% When relevant, we also include citations for each example, indicating the
% sources. The bibtex for each citation is included at the end of the
% library file.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%
% NOTE: The PIETOOLS implementation of 2D systems is currently very
% limited. PDEs have to specified in a very particular format, that we do
% not outline here. Feel free to play around with the 2D PIEs, but please
% wait for a later release for more information and enhanced functionality.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 28/09/2021
% Update for nargin=0 case, DJ - 12/23/21
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(2,'\n Warning: The 2D PIE implementation is still under construction. Functionality and documentation will be limited.\n')

% NOTE: The currently defined PDE struct will be cleared to avoid conflict
pvar ss1 ss2 tt1 tt2

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

BATCH = 0;      % if nonzero, batch-based PDE is assigned as output number BATCH of this function
TERM = 0;       % if nonzero, term-based PDE is assigned as output number TERM of this function
GUI = 0;        % 1 to load GUI, 0 not to

% Suppress (some) warnings in PDE initialization
evalin('base','silent_initialize_pde = true;');
% Each example comes with its own executives, so initially disable all
evalin('base','stability = 0;');
evalin('base','stability_dual = 0;');
evalin('base','Hinf_gain = 0;')
evalin('base','Hinf_gain_dual = 0;')
evalin('base','Hinf_estimator = 0;')
evalin('base','Hinf_control = 0;')

params = {};

% Collect the inputs
if nargin==0 %<-- If there is no input, pause the script and let the user specify an example
    userinp = input('\n Select an example (1 through 5) to convert \n ---> ','s');
    varargin0 = split(userinp,[" ",","]);
    index = str2double(varargin0{1});
    if ~isnan(index) && (index>=1 && index<=10)
        varargin0{1} = index;
    else
        userinp = input('\n No existing example specified... Please input an integer value 1 through 5 to extract the example \n ---> ','s');
        index = str2double(varargin0{1});
        if isnan(index)
            error('Please specify the desired example as first argument when calling ''examples_PDE_library_PIETOOLS''')
        else
            varargin0 = [str2double(userinp);varargin];
        end
    end
    nargin0 = length(varargin0);
    fprintf(['\n ... Extracting ODE-PDE example ', num2str(varargin0{1}),' ...\n']);
else
    nargin0 = nargin;
    varargin0 = varargin;
end

if nargin0==1 %<-- If there is one input, this input MUST correspond to the desired PDE
    if isdouble(varargin0{1}) && varargin0{1}<=10
        index = varargin0{1};
        TERM = 1;
    elseif contains(varargin0{1},'batch','IgnoreCase',true)
        index = randi(5,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        BATCH = 1;
    elseif contains(varargin0{1},'term','IgnoreCase',true)
        index = randi(5,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        TERM = 1;
    elseif contains(varargin0{1},'gui','IgnoreCase',true)
        index = randi(5,1);
        disp(['Warning: No example has been selected, randomly selecting example ',num2str(index)]);
        GUI = 1;
    else
        error('First argument must correspond to an example listed in this file')
    end
elseif nargin0>=2
    if isdouble(varargin0{1}) && varargin0{1}<=10
        index = varargin0{1};
    else
        error('First argument must correspond to an example listed in this file')
    end
    pindx = 1;
    for j=2:nargin0
        if contains(varargin0{j},'batch','IgnoreCase',true) || (isdouble(varargin0{j}) && (varargin0{j}==1))
            if BATCH==0
                BATCH = TERM+1;
            end
        elseif contains(varargin0{j},'term','IgnoreCase',true) || (isdouble(varargin0{j}) && (varargin0{j}==2))
            if TERM==0
                TERM = BATCH+1;
            end
        elseif contains(varargin0{j},'gui','IgnoreCase',true) || (isdouble(varargin0{j}) && (varargin0{j}==3))
            GUI = 1;
        elseif ischar(varargin0{j})
            try eval(varargin0{j});      %<-- In this case we assume the input defines certain parameters
                if contains(varargin0{j},';')
                    params{pindx} = varargin0{j};
                else
                    params{pindx} = [varargin0{j},';'];
                end
                pindx = pindx+1;
            catch
                disp(['Warning: Argument ',num2str(j),' is not understood, and is therefore ignored']);
            end
        else
            disp(['Warning: Argument ',num2str(j),' is not understood, and is therefore ignored']);
        end
    end
end
if BATCH==1 || GUI==1
    disp('Warning: Batch input format and GUI are not available for 2D systems. Reverting to terms input format...');
    TERM=1;
end
npars = length(params);

BATCH = 0;   GUI = 0;   TERM = 1;


if index==1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D PDEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %==================================================================
% % % 2D transport equation
% % % PDE         u_{t}  = cx*u_{x} + cy*u_{y}
% % % With BCs    u(x=0) = 0;   u(y=0) = 0;

 %%% Executive Function:
 evalin('base','stability = 1;');

 cx = 1; cy = 2; 
 n11 = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,1;0,1];
 PDE_t.n.n_pde = [0,0,0; 0,n11,0; 0,0,0];

 % PDE: u_{t} = cx * u_{x}
 PDE_t.PDE.A{1}.Lstate = [1,1];
 PDE_t.PDE.A{1}.Rstate = [1,1];
 PDE_t.PDE.A{1}.coeff = cx*eye(n11);
 PDE_t.PDE.A{1}.D = [1,0];

 % PDE: u_{t} = ... + cy * u_{y}
 PDE_t.PDE.A{2}.Lstate = [1,1];
 PDE_t.PDE.A{2}.Rstate = [1,1];
 PDE_t.PDE.A{2}.coeff = cy*eye(n11);
 PDE_t.PDE.A{2}.D = [0,1];

 BCindx = [ones(2*n11,1),zeros(2*n11,1);zeros(2*n11,2)];

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,n11,0,PDE_t.dom);



elseif index==2
% % % % %==================================================================
% % % 2D reaction-diffusion equation
% % % PDE         u_{t}  = lam*u + cx*u_{xx} + cy*u_{yy}
% % % With BCs    u(x=0) = 0;   u(y=0) = 0;
% % %             u(x=1) = 0;   u(y=1) = 0;
%
% % Stable when lam <= 2*pi^2 (Holmes, 1994 [1]).

 %%% Executive Function:
 evalin('base','stability = 1;');

 cx = 1; cy = 1; lam = 2*pi^2-1e-2;
 n22 = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,1;0,1];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];

 % PDE: u_{t} = lam * u
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = lam*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u_{t} = ... + cx * u_{xx}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = cx*eye(n22);
 PDE_t.PDE.A{2}.D = [2,0];

 % PDE: u_{t} = ... + cy * u_{yy}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = cy*eye(n22);
 PDE_t.PDE.A{3}.D = [0,2];

 BCindx = [ones(4*n22,1),zeros(4*n22,1)];

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);



elseif index==3
% % % % %==================================================================
% % % 2D parabolic equation
% % % PDE         u_{t}  = a*u + bx*u_{x} + by*u_{y} + cx*u_{xx} + cy*u_{yy}
% % % With BCs    u(x=0) = 0;       u(y=0) = 0;
% % %             u_{x}(x=0) = 0;   u_{y}(y=0) = 0;

 %%% Executive Function:
 evalin('base','stability = 1;');

 a = 1;  bx = 1;  by = 1;  cx = 1;  cy = 1;
 n22 = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,1;0,1];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];

 % PDE: u_{t} = a * u
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = a*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u_{t} = ... + bx * u_{x}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = bx*eye(n22);
 PDE_t.PDE.A{2}.D = [1,0];

 % PDE: u_{t} = ... + by * u_{y}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = by*eye(n22);
 PDE_t.PDE.A{3}.D = [0,1];

 % PDE: u_{t} = ... + cx * u_{xx}
 PDE_t.PDE.A{4}.Lstate = [2,2];
 PDE_t.PDE.A{4}.Rstate = [2,2];
 PDE_t.PDE.A{4}.coeff = cx*eye(n22);
 PDE_t.PDE.A{4}.D = [2,0];

 % PDE: u_{t} = ... + cy * u_{yy}
 PDE_t.PDE.A{5}.Lstate = [2,2];
 PDE_t.PDE.A{5}.Rstate = [2,2];
 PDE_t.PDE.A{5}.coeff = cy*eye(n22);
 PDE_t.PDE.A{5}.D = [0,2];

 BCindx = [ones(2*n22,2);zeros(2*n22,2)];

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);



elseif index==4
% % % % %==================================================================
% % % 2D wave equation
% % % PDE         u_{tt}  = cx*u_{xx} + cy*u_{yy};
% % % With BCs    u(x=0) = 0;       u(y=0) = 0;
% % %             u(x=1) = 0;       u(y=1) = 0;
%
% % We use states u1 = u, u2 = u_{t}.
% % Then u1(x=0) = u2(x=0) = u1(x=1) = u2(x=1) = 0,
% %      u1(y=0) = u2(y=0) = u1(y=1) = u2(y=1) = 0.

 %%% Executive Function:
 evalin('base','stability = 1;');

 cx = 1;  cy = 1;
 n22 = 2;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,1;0,1];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];

 % PDE: u1_{t} = u2
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = [0,1;0,0];
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u2_{t} = cx * u1_{xx}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = [0,0;cx,0];
 PDE_t.PDE.A{2}.D = [2,0];

 % PDE: u2_{t} = ... + cy * u1_{yy}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = [0,0;cy,0];
 PDE_t.PDE.A{3}.D = [0,2];

 BCindx = [ones(4*n22,1),zeros(4*n22,1)];

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);


 
elseif index==5
% % % % %==================================================================
% % % 2D heat equation coupled to ODE
% % % ODE         X_{t}  = (A+BK)*X + B*u(x=0,y=0)
% % % PDE         u_{t}  = cx*u_{xx} + cy*u_{yy}
% % % With BCs    u_{x}(x=0) = 0;   u_{y}(y=0) = 0;
% % %             u(x=1) = 0;       u(y=1) = 0;

 %%% Executive Function:
 evalin('base','stability = 1;');

 cx = 1;     cy = 1;     k = -2;
 nx = 1; n22 = 1;

if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end
 
 A = eye(nx);    B = eye(nx,n22);  K = k*eye(n22,nx);

 PDE_t.dom = [0,1;0,1];
 PDE_t.n.nx = 1;     PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];

 % ODE: X_{t} = (A+B*K) * X
 PDE_t.ODE.A = A+B*K;

 % ODE: X_{t} = ... + B*r
 PDE_t.ODE.Bxr = B;

 % PDE: r = u(0,0)
 PDE_t.PDE.Drb{1}.coeff = eye(n22); PDE_t.PDE.Drb{1}.Rstate = [2,2]; 
 PDE_t.PDE.Drb{1}.delta = [0,0]; PDE_t.PDE.Drb{1}.D = [0,0];  
 
 % PDE: u_{t} = cx * u_{xx}
 PDE_t.PDE.A{1}.Lstate = [2,2];  PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = cx*eye(n22);
 PDE_t.PDE.A{1}.D = [2,0];

 % PDE: u_{t} = ... + cy * u_{yy}
 PDE_t.PDE.A{2}.Lstate = [2,2];  PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = cy*eye(n22);
 PDE_t.PDE.A{2}.D = [0,2];

 % BCs: 0 = Ebb*u_{bf}
 BCindx = [zeros(2*n22,1),ones(2*n22,1);ones(2*n22,1),zeros(2*n22,1)];
 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);
 


elseif index==6
% % % % %==================================================================
% % % 2D reaction-diffusion equation with disturbance
% % % PDE         X_{t}  = lam*X + cx*X_{xx} + cy*X_{yy} + Bw*w(t)
% % % Output      z(t)   = int_dom Cz*X(t,x,y) dxdy
% % % With BCs    u(x=0) = 0;   u(y=0) = 0;
% % %             u(x=1) = 0;   u(y=1) = 0;
%
% % Stable for w(t)=0 when lam <= 2*pi^2 (Holmes, 1994 [1]).

 %%% Executive Function:
 evalin('base','Hinf_gain = 1;');

 n22 = 1;   nw = 1;     nz = 1;
 Lx = 1;    Ly = 1;
 cx = 1;    cy = 1;     lam = (pi/Lx)^2 + (pi/Ly)^2;
 Bw = 1*eye(n22,nw);    Cz = 1*eye(nz,n22);
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,Lx;0,Ly];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];
 PDE_t.n.nw = nw;   PDE_t.n.nz = nz;

 % ODE: v(t) = Bw*w(t)
 PDE_t.ODE.Dvw = Bw;
 
 % ODE: z(t) = r(t)
 PDE_t.ODE.Dzr = eye(nz);

 % PDE: u_{t} = lam * X
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = lam*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u_{t} = ... + cx * X_{xx}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = cx*eye(n22);
 PDE_t.PDE.A{2}.D = [2,0];

 % PDE: u_{t} = ... + cy * X_{yy}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = cy*eye(n22);
 PDE_t.PDE.A{3}.D = [0,2];
 
 % PDE: X_{t} = ... + v(t) = ... + Bw*w(t)
 PDE_t.PDE.Bpv = eye(n22);
 
 % PDE: z(t) = r(t) = int_{dom} Cz(x,y)*X(t,x,y) dx dy
 PDE_t.PDE.Crp{1}.Rstate = [2,2];
 PDE_t.PDE.Crp{1}.coeff = Cz;
 PDE_t.PDE.Crp{1}.D = [0,0];

 BCindx = [ones(4*n22,1), zeros(4*n22,1)]; % Dirichlet conditions on all boundaries

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);
 
 
 
elseif index==7
% % % % %==================================================================
% % % Custom parabolic PDE with all Dirichlet BCs
% % % PDE         X_{t}  = al(x,y)*(X_{xx} + X_{yy}) 
% % %                       + bet(x,y)*(X_{x}+X_{y}) + gam(x,y)*X 
% % %                           + Bw(x,y)*w(t)
% % % Output      z(t)   = int_dom Cz(x,y)*X(t,x,y) dxdy
% % % With BCs    u(x=0) = 0;   u(y=0) = 0;
% % %             u(x=1) = 0;   u(y=1) = 0;
%

 %%% Executive Function:
 evalin('base','Hinf_gain = 1;');

 Lx = 1;    Ly = 1;
 al = 1;    bet = 0;    gam = (pi/Lx)^2 + (pi/Ly)^2;
 Bw = 1;    Cz = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 n22 = size(al,2);  nw = size(Bw,2);   nz = size(Cz,1);

 PDE_t.dom = [-Lx/2,Lx/2;-Ly/2,Ly/2];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];
 PDE_t.n.nw = nw;   PDE_t.n.nz = nz;

 % ODE: v(t) = w(t)
 PDE_t.ODE.Dvw = eye(nw);
 
 % ODE: z(t) = r(t)
 PDE_t.ODE.Dzr = eye(nz);

 % PDE: u_{t} = gam(x,y) * X
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = gam*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u_{t} = ... + bet(x,y) * X_{x}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = bet*eye(n22);
 PDE_t.PDE.A{2}.D = [1,0];

 % PDE: u_{t} = ... + bet(x,y) * X_{y}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = bet*eye(n22);
 PDE_t.PDE.A{3}.D = [0,1];
 
 % PDE: u_{t} = ... + al(x,y) * X_{xx}
 PDE_t.PDE.A{4}.Lstate = [2,2];
 PDE_t.PDE.A{4}.Rstate = [2,2];
 PDE_t.PDE.A{4}.coeff = al*eye(n22);
 PDE_t.PDE.A{4}.D = [2,0];

 % PDE: u_{t} = ... + al(x,y) * X_{yy}
 PDE_t.PDE.A{5}.Lstate = [2,2];
 PDE_t.PDE.A{5}.Rstate = [2,2];
 PDE_t.PDE.A{5}.coeff = al*eye(n22);
 PDE_t.PDE.A{5}.D = [0,2];
 
 % PDE: X_{t} = ... + v(t) = ... + Bw*w(t)
 PDE_t.PDE.Bpv = Bw;
 
 % PDE: z(t) = r(t) = int_{dom} Cz(x,y)*X(t,x,y) dx dy
 PDE_t.PDE.Crp{1}.Rstate = [2,2];
 PDE_t.PDE.Crp{1}.coeff = Cz;
 PDE_t.PDE.Crp{1}.D = [0,0];

 BCindx = [ones(4*n22,1), zeros(4*n22,1)]; % Dirichlet conditions on all boundaries

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);
 
 
 
elseif index==8
% % % % %==================================================================
% % % Custom parabolic PDE with Dirichlet and Neumann BCs
% % % PDE         X_{t}  = al(x,y)*(X_{xx} + X_{yy}) 
% % %                       + bet(x,y)*(X_{x}+X_{y}) + gam(x,y)*X 
% % %                           + Bw(x,y)*w(t)
% % % Output      z(t)   = int_dom Cz(x,y)*X(t,x,y) dxdy
% % % With BCs    u(x=0) = 0;       u(y=0) = 0;
% % %             u_{x}(x=1) = 0;   u_{y}(y=1) = 0;
%
% For al=1, bet=0, will be stable when gam <= 0.5*pi^2

 %%% Executive Function:
 evalin('base','Hinf_gain = 1;');

 Lx = 1;    Ly = 1;
 al = 1;    bet = 0;    gam = (pi/Lx)^2 + (pi/Ly)^2;
 Bw = 1;    Cz = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 n22 = size(al,2);  nw = size(Bw,2);   nz = size(Cz,1);

 PDE_t.dom = [-Lx/2,Lx/2;-Ly/2,Ly/2];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];
 PDE_t.n.nw = nw;   PDE_t.n.nz = nz;

 % ODE: v(t) = w(t)
 PDE_t.ODE.Dvw = eye(nw);
 
 % ODE: z(t) = r(t)
 PDE_t.ODE.Dzr = eye(nz);

 % PDE: u_{t} = gam(x,y) * X
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = gam*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u_{t} = ... + bet(x,y) * X_{x}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = bet*eye(n22);
 PDE_t.PDE.A{2}.D = [1,0];

 % PDE: u_{t} = ... + bet(x,y) * X_{y}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = bet*eye(n22);
 PDE_t.PDE.A{3}.D = [0,1];
 
 % PDE: u_{t} = ... + al(x,y) * X_{xx}
 PDE_t.PDE.A{4}.Lstate = [2,2];
 PDE_t.PDE.A{4}.Rstate = [2,2];
 PDE_t.PDE.A{4}.coeff = al*eye(n22);
 PDE_t.PDE.A{4}.D = [2,0];

 % PDE: u_{t} = ... + al(x,y) * X_{yy}
 PDE_t.PDE.A{5}.Lstate = [2,2];
 PDE_t.PDE.A{5}.Rstate = [2,2];
 PDE_t.PDE.A{5}.coeff = al*eye(n22);
 PDE_t.PDE.A{5}.D = [0,2];
 
 % PDE: X_{t} = ... + v(t) = ... + Bw*w(t)
 PDE_t.PDE.Bpv = Bw;
 
 % PDE: z(t) = r(t) = int_{dom} Cz(x,y)*X(t,x,y) dx dy
 PDE_t.PDE.Crp{1}.Rstate = [2,2];
 PDE_t.PDE.Crp{1}.coeff = Cz;
 PDE_t.PDE.Crp{1}.D = [0,0];

 BCindx = [ones(2*n22,1), zeros(2*n22,1); zeros(2*n22,1), ones(2*n22,1)]; % Dirichlet and Neumann conditions on opposing boundaries

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);
 
 
 
 elseif index==9
% % % % %==================================================================
% % % Randomly generated parabolic PDE
% % % PDE         X_{t}  = (5*x^2-15*x*y+13*y^2)*X_{xx} + X_{yy} 
% % %                           + (10*x-15*y)*X_{x} + (-15*x+26*y)*X_{y}
% % %                               + (-17*x^4+30*y+25*x^2+8*y^3+50*y^4)*X 
% % %                                   + 16*x*(x-1)*y*(y-1)*w(t)
% % % Output      z(t)   = int_dom Cz*X(t,x,y) dxdy
% % % With BCs    u(x=0) = 0;   u(y=0) = 0;
% % %             u(x=1) = 0;   u(y=1) = 0;
%
% % [2], Meyer 2015

 %%% Executive Function:
 evalin('base','Hinf_gain = 1;');

 n22 = 1;   nw = n22;   nz = n22;
 Lx = 1;    Ly = 1;
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,Lx;0,Ly];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];
 PDE_t.n.nw = nw;   PDE_t.n.nz = nz;

 % ODE: v(t) = w(t)
 PDE_t.ODE.Dvw = eye(n22,nw);
 
 % ODE: z(t) = r(t)
 PDE_t.ODE.Dzr = eye(nz);

 % PDE: u_{t} = (-17*x^4+30*y+25*x^2+8*y^3+50*y^4) * X
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = (-17*ss1^4+30*ss2+25*ss1^2+8*ss2^3+50*ss2^4)*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: u_{t} = ... + (10*x-15*y) * X_{x}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = (10*ss1-15*ss2)*eye(n22);
 PDE_t.PDE.A{2}.D = [1,0];
 
 % PDE: u_{t} = ... + (-15*x+26*y) * X_{y}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = (-15*ss1+26*ss2)*eye(n22);
 PDE_t.PDE.A{3}.D = [0,1];
 
 % PDE: u_{t} = ... + (5*x^2-15*x*y+13*y^2) * X_{xx}
 PDE_t.PDE.A{4}.Lstate = [2,2];
 PDE_t.PDE.A{4}.Rstate = [2,2];
 PDE_t.PDE.A{4}.coeff = (5*ss1^2-15*ss1*ss2+13*ss2^2)*eye(n22);
 PDE_t.PDE.A{4}.D = [2,0];

 % PDE: u_{t} = ... + X_{yy}
 PDE_t.PDE.A{5}.Lstate = [2,2];
 PDE_t.PDE.A{5}.Rstate = [2,2];
 PDE_t.PDE.A{5}.coeff = eye(n22);
 PDE_t.PDE.A{5}.D = [0,2];
 
 % PDE: X_{t} = ... + (2/Lx)^2*(2/Ly^2)*x*(x-Lx)*y*(y-Ly)*v(t)
 PDE_t.PDE.Bpv = (2/Lx)^2*(2/Ly)^2* ss1*(ss1-Lx)* ss2*(ss2-Ly)* eye(n22);
 
 % PDE: z(t) = r(t) = int_{dom} Cz(x,y)*X(t,x,y) dx dy
 PDE_t.PDE.Crp{1}.Rstate = [2,2];
 PDE_t.PDE.Crp{1}.coeff = eye(nz,n22);
 PDE_t.PDE.Crp{1}.D = [0,0];

 BCindx = [ones(4*n22,1), zeros(4*n22,1)]; % Dirichlet conditions on all boundaries

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);
 
 
 
elseif index==10
% % % % %==================================================================
% % % Parabolic PDE with disturbance and actuator input
% % % PDE         X_{t}  = al*(X_{xx}+X_{yy}) + bet_x*X_{x} + bet_y*X_{y}
% % %                           + lam*X + Bw(x,y)*w(t) + Bu(x,y)*u(t)
% % % Output      z(t)   = int_{dom} Cz(x,y)*X(t,x,y) dxdy
% % % With BCs    X(x=0) = 0;   X(y=0) = 0;
% % %             X(x=1) = 0;   X(y=1) = 0;
%
% % [3], Demetriou 2019

 %%% Executive Function:
 evalin('base','Hinf_gain = 1;');

 Lx = 1;    Ly = 1;
 n22 = 1;   nw = 1;     nu = 1;     nz = 1;
 al = 0.1;  betx = 0;   bety = 0;       lam = (pi/Lx)^2 + (pi/Ly)^2;
 Bw = 1*eye(n22,nw);    Bu = 1*eye(n22,nu);     Cz = 1*eye(nz,n22);
 
if npars~=0
%%% Specify potential parameters
 for j=1:npars
     eval(params{j});
 end
end

 PDE_t.dom = [0,Lx;0,Ly];
 PDE_t.n.n_pde = [0,0,0; 0,0,0; 0,0,n22];
 PDE_t.n.nw = nw;   PDE_t.n.nu = nu;   PDE_t.n.nz = nz;

 % ODE: v(t) = Bw*w(t) + Bu*u(t)
 PDE_t.ODE.Dvw = Bw;
 PDE_t.ODE.Dvu = Bu;
 
 % ODE: z(t) = r(t)
 PDE_t.ODE.Dzr = eye(nz);

 % PDE: X_{t} = gam * X
 PDE_t.PDE.A{1}.Lstate = [2,2];
 PDE_t.PDE.A{1}.Rstate = [2,2];
 PDE_t.PDE.A{1}.coeff = lam*eye(n22);
 PDE_t.PDE.A{1}.D = [0,0];

 % PDE: X_{t} = ... + betx * X_{x}
 PDE_t.PDE.A{2}.Lstate = [2,2];
 PDE_t.PDE.A{2}.Rstate = [2,2];
 PDE_t.PDE.A{2}.coeff = betx*eye(n22);
 PDE_t.PDE.A{2}.D = [1,0];

 % PDE: X_{t} = ... + bety * X_{y}
 PDE_t.PDE.A{3}.Lstate = [2,2];
 PDE_t.PDE.A{3}.Rstate = [2,2];
 PDE_t.PDE.A{3}.coeff = bety*eye(n22);
 PDE_t.PDE.A{3}.D = [0,1];

 % PDE: X_{t} = ... + al * X_{xx}
 PDE_t.PDE.A{4}.Lstate = [2,2];
 PDE_t.PDE.A{4}.Rstate = [2,2];
 PDE_t.PDE.A{4}.coeff = al*eye(n22);
 PDE_t.PDE.A{4}.D = [2,0];

 % PDE: X_{t} = ... + al * X_{yy}
 PDE_t.PDE.A{5}.Lstate = [2,2];
 PDE_t.PDE.A{5}.Rstate = [2,2];
 PDE_t.PDE.A{5}.coeff = al*eye(n22);
 PDE_t.PDE.A{5}.D = [0,2];
 
 % PDE: X_{t} = ... + v(t) = ... + Bu*u(t) + Bw*w(t)
 PDE_t.PDE.Bpv = eye(n22);
 
 % PDE: z(t) = r(t) = int_{dom} Cz(x,y)*X(t,x,y) dx dy
 PDE_t.PDE.Crp{1}.Rstate = [2,2];
 PDE_t.PDE.Crp{1}.coeff = Cz;
 PDE_t.PDE.Crp{1}.D = [0,0];

 BCindx = [ones(4*n22,1),zeros(2*n22,1)]; % Dirichlet conditions on all boundaries

 PDE_t.BC.Ebb = DN_BCs2opvar2d(BCindx,0,n22,PDE_t.dom);

end

% Check if the number of desired outputs is reasonable
if nargout>0 && BATCH==0 && TERM==0
    error('No PDE structure has been produced. Use "Get PDE Object" in the GUI or specify ''batch'' or ''terms'' to obtain corresponding structure')
elseif nargout==2 && max(BATCH,TERM)==1
    error('To obtain PDE struct in both batch and term-based format, submit both ''batch'' and ''terms'' as arguments');
elseif nargout>2
    error('At most 2 outputs can be produced: A PDE struct in batch format, and a PDE struct in term-based format');
end

% Specify the outputs of the function
if BATCH~=0
    varargout{BATCH} = PDE_b;
end
if TERM~=0
    varargout{TERM} = PDE_t;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % [1] - 
% E. E. Holmes, M. Lewis, J. Banks, and R. R. Veit, “Partial differential
% equations in ecology: Spatial interactions and population dynamics,”
% Ecology, vol. 75, pp. 17–29, 1994.
%
% % [2] -
% @inproceedings{meyer2015stability,
%   title={Stability analysis of parabolic linear PDEs with two spatial dimensions using Lyapunov method and SOS},
%   author={Meyer, Evgeny and Peet, Matthew M},
%   booktitle={2015 54th IEEE Conference on Decision and Control (CDC)},
%   pages={1884--1890},
%   year={2015},
%   organization={IEEE}
% }
%
% % [3] - 
% @inproceedings{demetriou2019feedback,
% title={Feedback kernel approximations and sensor selection for controlled 2D parabolic PDEs using computational geometry methods},
% author={Demetriou, Michael A and Hu, Weiwei},
% booktitle={2019 IEEE 58th Conference on Decision and Control (CDC)},
% pages={2144--2150},
% year={2019},
% organization={IEEE}
% }