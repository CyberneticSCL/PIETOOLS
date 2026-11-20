%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_options_check.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [opts, uinput]=PIESIM_options_check(varinput)
% Check if option and uinput fields are supplied by the user
% Initialze options to default options if not
% Initialize uinput as empty structure if not. The fields of uinput are not
% checked here.
%
% An argument after the first is opts if it has one of the option fields,
% uinput if it has one of the fields ic, w, exact. Option defaults: N = 8,
% tf = 1, intScheme = 1, Norder = 2, dt = 0.01, plot = 'no', ploteig =
% 'no', ifexact = false; intScheme outside {1,2} and Norder outside
% {1,...,4} are reset to their defaults, tf < dt (tf ~= 0) to tf = 0.
%
% NOTE: No other variables are checked here.
% Inputs:
% varinput - variable number of arguments (between 1 and 4)
% Required:
% 1) varinput(1): data structure of the proglem: PDE, DDE or PIE
% PIE structure of the problem specifies PI operators, T,Tu,Tw, A, Bi, Ci, Dij as fields
% Optional:
% 2) varinput(2): opts - options for simulation parameters. If empty or incomplete, will be
% set to default values
% 3) varinput(3): uinput - user-defined boundary inputs, forcing and initial
% conditions. If empty or incomplete, will be set to default values
%
% Outputs: 
% 1) opts - user-defined option fields (if defined), default option fields (if undefined). 
%    Also Nsteps = floor(tf/dt) (0 if dt = 0) and type: 'PIE', 'DDE',
%    'PDE_b' (legacy batch format), 'PDE_t' (legacy terms format) or 'PDE'
% 2) uinput - user's defined structure (if defined), empty structure (if
% undefiend)
% All properly defined variables are uchanged.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIESIM_options_check
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu,
% or Y. Peet at ypeet@asu.edu
% Initial coding YP  - 12_22_2021
% Added a clause -  "If piesize is not found, default to last argument" -
% 6/16/2022
% Update to account for new terms format - DJ, 10/10/2022.
% MMP, 11/19/2026: Release comment pass (comment-only).


% Number of arguments passed
nargin=length(varinput);

% required fields for options and uinputs
fields_opts = {'N','tf','intScheme','Norder','dt','plot','ploteig','ifexact'};
default_opts = {8, 1, 1, 2, 0.01,'no','no',false};
fields_uinput = {'ic','w','exact'};

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
if (opts.dt~=0)
opts.Nsteps=floor(opts.tf/opts.dt);
else
opts.Nsteps=0;
end
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

% Check if time step and final time are defined correctly

if (opts.tf<opts.dt & opts.tf~=0)
    disp('Warning: final time is less than the time step. Defaulting final time to zero');
    opts.tf=0;
end


%-------------------------------------
% Determine the type of the problem
%-------------------------------------

structure = varinput{1};

% PIE: 'pie_struct' or field T; DDE: field tau; legacy batch PDE: fields
% no, n0, n1 or n2; legacy terms PDE: fields n, ODE or PDE; PDE:
% 'pde_struct', 'sys', or fields BC or X
if isa(structure,'pie_struct') || isfield(structure,'T')
    disp('Solving PIE problem');
    opts.type='PIE';
elseif isfield(structure,'tau')
    disp('Solving DDE problem');
    opts.type='DDE';
elseif isfield(structure,'no') || isfield(structure,'n0') || isfield(structure,'n1') || isfield(structure,'n2')
    disp('Solving PDE problem in legacy batch format');
    opts.type='PDE_b';
elseif isfield(structure,'n') || isfield(structure,'ODE')  || isfield(structure,'PDE')
    disp('Solving PDE problem in legacy terms format');
    opts.type='PDE_t';
elseif isa(structure,'pde_struct') || isa(structure,'sys')||isfield(structure,'BC') || isfield(structure,'X')
    disp('Solving PDE problem');
    opts.type='PDE';
else
    error('Data structure must be specified')
end

    
end
