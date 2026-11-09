function varargout = examples_PDE_library_PIESOS(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% examples_PDE_library_PIESOS.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This library file contains the definition of some nonlinear PDE systems
% drawn from the literature. To use one of the examples, call this function
% with as FIRST input the index of the example you wish to use. Please
% check the table of contents below to see which index corresponds to which
% problem. The output of the function is the corresponding 'pde_struct'
% object describing the PDE.
%
% Inputs:
%     PDE = examples_PDE_library_PIESOS(example,option1,option2,...)
%       - example - integer (1-10) corresponding to the selected example
%       - option# - character array 'name=value', which overrides the
%            default value of the parameter 'name' of the selected example
%            by the value 'value'. For example, 'r=5' sets r equal to 5.
%            The value can be any Matlab expression, e.g. 'b=2' or
%            'r=pi^2-0.1'. The parameters that can be overridden for each
%            example are listed in the table of contents below.
%            Alternatively, option# can be 'light', 'default' or 'heavy',
%            specifying the settings PIESOS_settings(option#) with which to
%            run the stability test associated to the example (see below).
%            In that case, the test is run without asking for permission.
%
% Nonlinear PDEs can only be declared in the command-line (pde_var) input
% format; the batch, term and GUI input formats are not supported.
%
% The output PDE can be converted to a PIE, and analyzed using PIESOS, e.g.
%     PDE = examples_PDE_library_PIESOS(2,'r=5');
%     PIE = convert(PDE);
%
% Most examples come with a stability test that can be verified for the
% default values of the parameters, using the PIESOS executives. This test
% is specified by
%       use_exec - 'global' for PIESOS_stability_global, or 'local' for
%            PIESOS_stability_local ('' if no test could be verified);
%       lambda   - the decay rate at which to verify stability;
%       LF_form  - the form of the Lyapunov functional, 'pie' or 'pde';
%       rad, alpha - (local test only) the radius of the ball, and the
%            weights of the Sobolev norm defining the ball;
%       settings - the settings of the executive, the default settings
%            PIESOS_settings() unless specified otherwise by the example;
% After extracting the PDE, the library asks whether you wish to run this
% test, and with what settings: 'light', 'default' or 'heavy', or, if you
% just press enter, those specified by the example. If you change the
% parameters of the example, the test may no longer be feasible.
%
% If you wish to include a new example in our library, please send it to us
% and we will include it in the next release. Please also include the
% citation information, if available.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 11/08/2026: Initial coding;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process the user inputs
[index,params,sttngs] = process_inputs(varargin,nargin);

pvar t s s1 s2 s3

% Stability test to run for the example: executive, ball radius and norm
% weights (local test), decay rate, form of the Lyapunov functional, and
% settings of the executive
use_exec = '';      rad = [];       alpha = 1;
lambda = 0;         LF_form = 'pie';
settings = PIESOS_settings();

fprintf([' --- Extracting PDE example ', num2str(index),' ---\n']);
switch index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Table of Contents
%
%__________________________________________________________________________
%   Problem                                      | Parameters
%==========================================================================
%  Global stability tests
%   1. Nonlinear delay ODE                       | alp, bet, tau
%   2. Burgers' equation with reaction term      | r, a, b
%   3. Reaction-diffusion eq., cubic damping     | r, lam, a, b
%   4. Korteweg-de Vries equation                | a, b
%   5. Split Burgers' equation, reaction term    | r, a, b
%   6. Klein-Gordon eq., nonlinear damping       | a, b
%   7. Kuramoto-Sivashinsky equation             | r, a, b
%   8. 2D Burgers' equation, reaction term       | r, a, b
%  Local stability tests
%   9. Fisher equation                           | alp, bet, a, b
%   10. Nonlinear delay ODE, delay-dependent     | a, tau
%--------------------------------------------------------------------------
    case 1
%% Nonlinear delay ODE
% 	ODE:    x_{t} = -alp*x(t) + bet*x(t-tau) - x(t)^3
% The zero solution is stable for all tau if alp >= |bet|, and unstable for
% all tau if alp < bet, see [1].

alp = 2;                % decay rate
bet = 1;                % delayed gain
tau = 1;                % delay
[alp,bet,tau] = set_params({'alp','bet','tau'},{alp,bet,tau},params);

x = pde_var('state',1,[],[]);
PDE = diff(x,t)==-alp*x+bet*subs(x,t,t-tau)-x^3;

% Stability test: global stability, verified with default settings
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pie';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
    case 2
%% Burgers' equation with reaction term
% 	PDE:    x_{t} = x_{ss} + r*x - x*x_{s}                     s in [a,b]
%   BCs:    x(s=a) = x(s=b) = 0;
% The zero solution is globally stable iff r <= pi^2/(b-a)^2, the 
% instability threshold of the reaction-diffusion equation. Global 
% stability can be verified for r <= pi^2-0.001, for [a,b]=[0,1].

a = 0;                  % left end of the interval of s
b = 1;                  % right end of the interval of s
r = pi^2 - 0.001;       % reaction parameter
[r,a,b] = set_params({'r','a','b'},{r,a,b},params);

x = pde_var(s,[a,b]);
PDE = [diff(x,t)==diff(x,s,2)+r*x-x*diff(x,s);
        subs(x,s,a)==0;
        subs(x,s,b)==0];

% Stability test: global stability, verified with default settings
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pie';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
    case 3
%% Reaction-diffusion equation with cubic damping
% 	PDE:    x_{t} = x_{ss} + r*x - x^3                         s in [a,b]
%   BCs:    x(s=a) = x(s=b) = 0;
% The zero solution is globally exponentially stable at decay rate lam, for 
% any lam <= pi^2/(b-a)^2 - r, see [2]. For r=5 and [a,b]=[0,1], this
% rate is 4.8696, and global stability can be verified for lam=4.7.

a = 0;                  % left end of the interval of s
b = 1;                  % right end of the interval of s
r = 5;                  % reaction parameter
lam = 4.7;              % decay rate at which to verify stability
[r,lam,a,b] = set_params({'r','lam','a','b'},{r,lam,a,b},params);

x = pde_var(s,[a,b]);
PDE = [diff(x,t)==diff(x,s,2)+r*x-x^3;
        subs(x,s,a)==0;
        subs(x,s,b)==0];

% Stability test: global exponential stability at decay rate lam, verified
% with default settings
use_exec = 'global';    % executive to use
lambda = lam;           % decay rate
LF_form = 'pde';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
    case 4
%% Korteweg--de Vries equation
% 	PDE:    x_{t} = -x_{sss} - x*x_{s}                         s in [a,b]
%   BCs:    x(s=a) = x(s=b) = 0;    x_{s}(s=b) = 0;
% The L2 norm of x does not increase, as d/dt ||x||^2 = -x_{s}(a)^2. Global
% stability can be verified for [a,b]=[0,1].

a = 0;                  % left end of the interval of s
b = 1;                  % right end of the interval of s
[a,b] = set_params({'a','b'},{a,b},params);

x = pde_var(s,[a,b]);
PDE = [diff(x,t)==-diff(x,s,3)-x*diff(x,s);
        subs(x,s,a)==0;
        subs(x,s,b)==0;
        subs(diff(x,s),s,b)==0];

% Stability test: global stability, verified with default settings
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pie';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
    case 5
%% Split Burgers' equation with reaction term
% 	PDE:    x1_{t} = x1_{ss} + r*x1 - x1*x1_{s}
%           x2_{t} = x2_{ss} + r*x2 - x2*x2_{s}         s in [a,(a+b)/2]
%   BCs:    x1(s=a) = 0;        x2(s=(a+b)/2) = 0;
%           x1(s=(a+b)/2) = x2(s=a);
%           x1_{s}(s=(a+b)/2) = x2_{s}(s=a);
% The zero solution is globally stable iff r <= pi^2/(b-a)^2. Global 
% stability can be verified for r=pi^2-0.01 and [a,b]=[0,1].

a = 0;                  % left end of the interval of s
b = 1;                  % right end of the interval of s
r = pi^2 - 0.01;        % reaction parameter
[r,a,b] = set_params({'r','a','b'},{r,a,b},params);

c = (a+b)/2;            % midpoint, where the interval is split
x1 = pde_var(s,[a,c]);
x2 = pde_var(s,[a,c]);
PDE = [diff(x1,t)==diff(x1,s,2)+r*x1-x1*diff(x1,s);
        diff(x2,t)==diff(x2,s,2)+r*x2-x2*diff(x2,s);
        subs(x1,s,a)==0;
        subs(x2,s,c)==0;
        subs(x1,s,c)==subs(x2,s,a);
        subs(diff(x1,s),s,c)==subs(diff(x2,s),s,a)];

% Modified settings for global stability test
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pde';        % form of the Lyapunov functional
settings.LF.mon_deg = 2;                % degree of the kernels of V
settings.LF.exclude = [0;0;0];          % allow integral terms in V
settings.LF.sep = true;                 % separable kernels in V
settings.sos.mon_deg = [NaN,NaN,6];     % degree of the derivative constraint
settings.sos.psatz = [0,0,1];           % psatz multipliers in the constraints

%--------------------------------------------------------------------------
    case 6
%% Klein-Gordon equation with nonlinear damping
% 	PDE:    u_{tt} = u_{ss} - u_{t}*u_{s}^2                    s in [a,b]
%   BCs:    u(s=a) = 0;     u_{s}(s=b) = 0;
% Written in first order form, using x1 = u_{t} and x2 = u_{s}:
% 	PDE:    x1_{t} = x2_{s} - x1*x2^2;
%           x2_{t} = x1_{s};
%   BCs:    x1(s=a) = 0;    x2(s=b) = 0;
% Global stability can be verified for [a,b]=[0,1].

a = 0;                  % left end of the interval of s
b = 1;                  % right end of the interval of s
[a,b] = set_params({'a','b'},{a,b},params);

x1 = pde_var(s,[a,b]);
x2 = pde_var(s,[a,b]);
PDE = [diff(x1,t)==diff(x2,s)-x1*x2^2;
        diff(x2,t)==diff(x1,s);
        subs(x1,s,a)==0;
        subs(x2,s,b)==0];

% Stability test: global stability, verified with default settings.
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pde';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
    case 7
%% Kuramoto-Sivashinsky equation
% 	PDE:    x_{t} = -x_{ssss} - x_{ss} - x*x_{s}       s in [a,b]
%   BCs:    x(s=a) = x(s=b) = 0;    x_{s}(s=a) = x_{s}(s=b) = 0;
% The zero solution is globally exponentially stable iff b-a < 2*pi, the 
% buckling length of the linearization (see [3]). Global stability can be 
% verified for [a,b]=[0,5.5], using a Lyapunov functional with constant 
% kernel.

a = 0;                  % left end of the interval of s
b = 5.5;                % right end of the interval of s
[a,b] = set_params({'a','b'},{a,b},params);

x = pde_var(s,[a,b]);
PDE = [diff(x,t)==-diff(x,s,4)-diff(x,s,2)-x*diff(x,s);
        subs(x,s,a)==0;    subs(diff(x,s),s,a)==0;
        subs(x,s,b)==0;    subs(diff(x,s),s,b)==0];

% Modified settings to reduce computation time
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pde';        % form of the Lyapunov functional
settings.LF.mon_deg = 0;                % constant kernel of V
settings.sos.psatz = [0,0,0];           % no psatz multipliers

%--------------------------------------------------------------------------
    case 8
%% 2D Burgers' equation with reaction term
% 	PDE:    x_{t} = x_{s1s1} + x_{s2s2} + r*x - x*(x_{s1} + x_{s2});
%   BCs:    x = 0 on the boundary of [a,b]^2;               s1,s2 in [a,b]
% The zero solution is globally stable iff r <= 2*pi^2/(b-a)^2. 
% Global stability can be verified for r=15 and [a,b]=[0,1].

a = 0;                  % left end of the intervals of s1 and s2
b = 1;                  % right end of the intervals of s1 and s2
r = 15;                 % reaction parameter
[r,a,b] = set_params({'r','a','b'},{r,a,b},params);

x = pde_var([s1,s2],[a,b;a,b]);
PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+r*x-x*diff(x,s1)-x*diff(x,s2);
        subs(x,s1,a)==0;    subs(x,s1,b)==0;
        subs(x,s2,a)==0;    subs(x,s2,b)==0];

% Modified settings to reduce computation time
use_exec = 'global';    % executive to use
lambda = 0;             % decay rate
LF_form = 'pde';        % form of the Lyapunov functional
settings.LF.mon_deg = 1;                % degree of the kernels of V
settings.sos.mon_deg = [1,1,3];         % degrees of the SOS constraints
settings.sos.psatz = [0,0,0];           % no psatz multipliers

%--------------------------------------------------------------------------
    case 9
%% Fisher equation
% 	PDE:    x_{t} = x_{ss} + alp*x - bet*x^2                   s in [a,b]
%   BCs:    x(s=a) = x(s=b) = 0;
% The zero solution is locally stable iff alp <= pi^2/(b-a)^2 (see [4]). 
% Local stability can be verified on the L2 ball of radius 2, for alp=5, 
% bet=-1, and [a,b]=[0,1].

a = 0;                  % left end of the interval of s
b = 1;                  % right end of the interval of s
alp = 5;                % linear reaction parameter
bet = -1;               % quadratic reaction parameter
[alp,bet,a,b] = set_params({'alp','bet','a','b'},{alp,bet,a,b},params);

x = pde_var(s,[a,b]);
PDE = [diff(x,t)==diff(x,s,2)+alp*x-bet*x^2;
        subs(x,s,a)==0;
        subs(x,s,b)==0];

% Stability test: local stability on the L2 ball of radius 2, verified
% with default settings
use_exec = 'local';     % executive to use
rad = 2;                % radius of the ball
alpha = 1;              % weights of the norm: 1 for the L2 norm
lambda = 0;             % decay rate
LF_form = 'pie';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
    case 10
%% Nonlinear delay ODE, delay-dependent stability
% 	ODE:    x_{t} = -a*x(t-tau) - x(t)^3
% The zero solution is locally stable if a*tau < pi/2, and unstable if 
% a*tau > pi/2. For tau=1, local stability can be verified on the L2 ball 
% of radius 0.5 (in x(t) and its history) up to a=1.55. 

a = 1.5;                % gain
tau = 1;                % delay
[a,tau] = set_params({'a','tau'},{a,tau},params);

x = pde_var('state',1,[],[]);
PDE = diff(x,t)==-a*subs(x,t,t-tau)-x^3;

% Stability test: local stability on the L2 ball of radius 0.5, verified
% with default settings for tau=1
use_exec = 'local';     % executive to use
rad = 0.5;              % radius of the ball
alpha = 1;              % weights of the norm: 1 for the L2 norm
lambda = 0;             % decay rate
LF_form = 'pie';        % form of the Lyapunov functional

%--------------------------------------------------------------------------
end

% % % Processing of the outputs.
% Initialize the PDE
PDE = initialize(PDE);
display(PDE);
varargout{1} = PDE;

% Check if the user wants to run the stability test associated to the example
if isempty(use_exec)
    fprintf('\n No stability test could be verified for this example.\n');
    return
end
if ~isempty(fieldnames(params))
    fprintf('\n The stability test was verified for the default parameter values only.\n');
end
if isempty(sttngs)
    % No settings specified by the user: ask whether to run the test, and
    % how hard the solver should work
    userinp = input('\n Would you like to run the stability test associated to this problem? (y/n) \n ---> ','s');
    userinp = strrep(userinp,'''','');    % Get rid of potential apostrophes
    if ~any(strcmpi(userinp,{'y','yes','n','no'}))
        userinp = input(' Please answer y (yes) or n (no). \n ---> ','s');
        userinp = strrep(userinp,'''','');
        if ~any(strcmpi(userinp,{'y','yes','n','no'}))
            error("Response not recognized; execution of the stability test aborted.")
        end
    end
    if any(strcmpi(userinp,{'n','no'}))
        return
    end
    fprintf('\n How hard should the solver work? \n');
    msg = ['   Please input ''light'', ''default'', or ''heavy'', or press enter to use the settings specified for this example \n ---> '];
    sttngs = input(msg,'s');
    sttngs = strrep(sttngs,'''','');
    if ~isempty(sttngs) && ~any(strcmpi(sttngs,{'light','default','heavy'}))
        fprintf(2,"\n Unknown settings specified; proceeding with the settings specified for this example.\n")
        sttngs = '';
    end
end
if ~isempty(sttngs)
    settings = PIESOS_settings(lower(sttngs));
    fprintf("\n Proceeding with '%s' settings.\n",lower(sttngs))
end

% Run the stability test with the specified values
PIE = convert(PDE);
if strcmp(use_exec,'global')
    PIESOS_stability_global(PIE,lambda,LF_form,settings);
else
    PIESOS_stability_local(PIE,rad,alpha,lambda,LF_form,settings);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [index,params,sttngs] = process_inputs(varargin0,nargin1)
% [INDEX,PARAMS,STTNGS] = PROCESS_INPUTS(VARARGIN0,NARGIN1) processes the
% user inputs of the examples_PDE_library_PIESOS function.
%
% INPUT
% - varargin0:  1 x nargin1 cell with the inputs of the library function.
%               The first is the index of the example, the others
%               character arrays 'name=value', each overriding the default
%               value of the parameter 'name' by the value of the expression
%               'value', or one of 'light', 'default', 'heavy', specifying
%               the settings of the stability test. If nargin1=0, the index
%               is requested from the user;
% - nargin1:    number of inputs of the library function;
%
% OUTPUT
% - index:      integer index of the example;
% - params:     struct with a field for each parameter specified by the
%               user, holding the specified value. Has no fields if no
%               parameters are specified. See 'set_params';
% - sttngs:     'light', 'default' or 'heavy' if the user specified the
%               settings of the stability test, and '' otherwise;

n_examples = 10;

params = struct();      % no fields if the user specifies no parameters
sttngs = '';

% Determine the index of the example, asking the user if none is specified
if nargin1==0
    userinp = input(['\n Select an example (1 through ',num2str(n_examples),') to extract \n ---> '],'s');
    index = str2double(userinp);
    varargin0 = {index};
    nargin1 = 1;
else
    index = varargin0{1};
end
if ~isnumeric(index) || ~isscalar(index) || index~=round(index) || index<1 || index>n_examples
    error(['First argument must be an integer 1 through ',num2str(n_examples),', the index of an example listed in this file.'])
end

% Process the settings and parameter specifications 'name=value'
for j=2:nargin1
    arg = varargin0{j};
    if isstring(arg) && isscalar(arg)
        arg = char(arg);
    end
    k = [];
    if ischar(arg)
        if any(strcmpi(strtrim(arg),{'light','default','heavy'}))
            sttngs = lower(strtrim(arg));
            continue
        end
        k = find(arg=='=',1);
    end
    if isempty(k)
        warning(['Argument ',num2str(j),' is not understood (parameters must be specified as ''name=value'', and settings as ''light'', ''default'' or ''heavy''), and is therefore ignored.']);
        continue
    end
    name = strtrim(arg(1:k-1));
    if ~isvarname(name)
        warning(['Argument ',num2str(j),' is not understood (''',name,''' is not a valid parameter name), and is therefore ignored.']);
        continue
    end
    try
        % The value may be any expression, e.g. pi^2-0.1 or [0,2]
        params.(name) = eval(regexprep(strtrim(arg(k+1:end)),';+$',''));
    catch
        warning(['Argument ',num2str(j),' is not understood (the value of ''',name,''' cannot be evaluated), and is therefore ignored.']);
    end
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function varargout = set_params(names,vals,params)
% VARARGOUT = SET_PARAMS(NAMES,VALS,PARAMS) overwrites the default values
% VALS of the parameters NAMES of an example by the values specified by the
% user, e.g. [r,a,b] = set_params({'r','a','b'},{r,a,b},params).
%
% INPUT
% - names:      1 x n cellstr with the names of the parameters of the example;
% - vals:       1 x n cell with the default values of the parameters;
% - params:     struct with a field for each parameter specified by the
%               user, as output by 'process_inputs';
%
% OUTPUT
% - varargout:  n outputs, the value params.(names{i}) if the user specified
%               parameter i, and the default value vals{i} otherwise.
%               A warning is thrown for each field of params that is not in
%               names, which is then ignored.

if ~iscellstr(names) || ~iscell(vals) || numel(names)~=numel(vals)
    error('Parameter names and default values must be specified as cell arrays of the same size.')
end
varargout = vals(:)';

% Overwrite the defaults by the values specified by the user
for i=1:numel(names)
    if isfield(params,names{i})
        varargout{i} = params.(names{i});
    end
end

% Warn for specified parameters that this example does not have
unknown = setdiff(fieldnames(params),names(:));
for i=1:numel(unknown)
    warning(['''',unknown{i},''' is not a parameter of this example, and is therefore ignored. Parameters: ',strjoin(names,', '),'.']);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% References %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% % [1] - Nonlinear delay ODE: delay-independent stability of the scalar
% %       delay equation
% @book{hale1993introduction,
%   title={Introduction to Functional Differential Equations},
%   author={Hale, Jack K. and Verduyn Lunel, Sjoerd M.},
%   series={Applied Mathematical Sciences},
%   volume={99},
%   year={1993},
%   publisher={Springer}
% }
%
% % [2] - Reaction-diffusion equation with cubic damping: equilibria
% %       bifurcating from zero at r = n^2*pi^2/(b-a)^2
% @article{chafee1974bifurcation,
%   title={A bifurcation problem for a nonlinear partial differential
%          equation of parabolic type},
%   author={Chafee, Nathaniel and Infante, Ettore F.},
%   journal={Applicable Analysis},
%   volume={4},
%   number={1},
%   pages={17--37},
%   year={1974}
% }
%
% % [3] - Kuramoto-Sivashinsky equation: stability iff lambda < 4*pi^2
% @article{liu2001stability,
%   title={Stability enhancement by boundary control in the
%          {K}uramoto--{S}ivashinsky equation},
%   author={Liu, Wei-Jiu and Krstic, Miroslav},
%   journal={Nonlinear Analysis: Theory, Methods \& Applications},
%   volume={43},
%   number={4},
%   pages={485--507},
%   year={2001}
% }
%
% % [4] - Fisher equation: critical patch size for absorbing boundaries
% @article{skellam1951random,
%   title={Random dispersal in theoretical populations},
%   author={Skellam, J. G.},
%   journal={Biometrika},
%   volume={38},
%   number={1-2},
%   pages={196--218},
%   year={1951}
% }
%
