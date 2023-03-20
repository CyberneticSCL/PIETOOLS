%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_options_check.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts, uinput]=PIESIM_options_check(varinput)
% Check if option and uinput fields are supplied by the user
% Initialze options to default options if not
% Initialize uinput as empty structure if not. Check and re-initialization of uinput fields will be
% done later in PIESIM_input_check

% NOTE: All other variables will be checked in PIETOOLS converter
% Inputs:
% varinput - variable number of arguments (between 1 and 4)
% Required:
% 1) varinput(1): data structure of the proglem: PDE, DDE or PIE
% PIE structure of the problem specifies PI operators, T,Tu,Tw, A, Bi, Ci, Dij as fields
% if varargin(1) is PDE or DDE, the rest of the inputs are optional
% if varargin(1) is PIE, the rest of the inputs are requires
% 2) varargin(2): opts - options for simulation parameters. If empty or incomplete, will be
% set to default values
% 3) varargin(3): uinput - user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values
% Not used for PDE/DDE, required for PIE
% 4) varargin(4): n_pde - number of states with increasing differentiability, for example
% [1,2,3] stands for (1) continuous state, (2) continuously differentiable,
% and (3) twice continuously differentiable states - only used it data structure is PIE 

% Outputs: 
% 1) opts - user-defined option fields (if defined), default option fields (if undefined). 
% 2) uinput - user's defined structure (if defined), empty structure (if
% undefiend)
% All properly defined variables are uchanged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 12_22_2021
% Added a clause -  "If piesize is not found, default to last argument" -
% 6/16/2022
% Update to account for new terms format - DJ, 10/10/2022.


nargin=length(varinput);

% required fields for options and uinputs
fields_opts = {'N','tf','intScheme','Norder','dt','plot','ploteig'};
default_opts = {8, 1, 1, 2, 0.01,'no','no'};
fields_uinput = {'ic','w','ifexact','exact'};

kopt=0;
kuopt=0;

if nargin==1
    opts = struct();
    uinput = struct();
else
    % Check if options field is passed
    for k=2:nargin
        optdefine=0;
        for i=1:length(fields_opts)
            if isfield(varinput{k},fields_opts{i})
                optdefine=1;
                kopt=k;
                break;
            end
        end
        if (optdefine==1)
            break
        end
    end
    if (optdefine==1)
        opts = varinput{kopt};
    else
        opts=struct();
    end
    
    % Check if uinput field is passed
    for k=2:nargin
        optdefine=0;
        for i=1:length(fields_uinput)
            if isfield(varinput{k},fields_uinput{i})
                optdefine=1;
                kuopt=k;
                break
            end
        end
        if (optdefine==1)
            break
        end
    end
    
    if (optdefine==1)
        uinput = varinput{kuopt};
    else
        uinput=struct();
    end
    
end

% check if all options are defined, if not define them
for i=1:length(fields_opts)
    if ~isfield(opts,fields_opts{i})
        opts.(fields_opts{i}) = default_opts{i};
        X = ['Warning: option  ',fields_opts{i},' is not defined. Setting to a default value of ', num2str(default_opts{i})];
        disp(X);
    end
end
opts.Nsteps=floor(opts.tf/opts.dt);
% Check if opts.intScheme is defined correctly
if (opts.intScheme~=1&opts.intScheme~=2)
    i=3;
    opts.(fields_opts{i}) = default_opts{i};
        X = ['Warning: option ',fields_opts{i},'  is out of bounds. Setting to a default value of ', num2str(default_opts{i})];
    disp(X);
end


% Check if opts.Norder is defined correctly
if (~ismember(opts.Norder,[1,2,3,4]))
    i=4;
    opts.(fields_opts{i}) = default_opts{i};
        X = ['Warning: option ',fields_opts{i},'  is out of bounds. Setting to a default value of ', num2str(default_opts{i})];
    disp(X);
end


%-------------------------------------
% Determine the type of the problem
%-------------------------------------

structure = varinput{1};

if isa(structure,'pie_struct') || isfield(structure,'T')
    disp('Solving PIE problem');
    opts.type='PIE';
    
    kargs = 2:nargin;
    ksize = kargs(find(kargs~=kopt));
    ksize = ksize(find(ksize~=kuopt));
    if (isempty(ksize))
        error("If input object type is PIE, then number of differentiable states should be specified, for example 'executive_PIESIM(PIE,opts,uinput,n_pde)'");
    else
        opts.piesize=varinput{ksize};
        % If piesize is not found, default to last argument
        if (isempty(opts.piesize))
            opts.piesize=varinput{nargin};
        end
    end
elseif isfield(structure,'tau')
    disp('Solving DDE problem');
    opts.type='DDE';
elseif isa(structure,'pde_struct') || isa(structure,'sys')
    disp('Solving PDE problem');
    opts.type='PDE';
elseif isfield(structure,'nx') || isfield(structure,'n0') || isfield(structure,'n1') || isfield(structure,'n2')
    disp('Solving PDE problem in legacy batch format');
    opts.type='PDE_b';
elseif isfield(structure,'n') || isfield(structure,'ODE')  || isfield(structure,'PDE')
    disp('Solving PDE problem in legacy terms format');
    opts.type='PDE_t';
else
    error('Data structure must be specified')
end
    
end