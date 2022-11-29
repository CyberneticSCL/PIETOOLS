%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_AUTO_EXECUTE.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A script to aid users in analyzing and controlling PIEs using PIETOOLS
% If a desired executive has been specified (e.g. 'Hinf_gain=1'), and a
% struct of settings is available in the workspace, this script will run
% the desired executive with the desired settings.
% If no executive has been specified, the script will ask the user which
% executive they want to run.
% If no settings have been specified, the script will ask the user what
% setting they want to use.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2021  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 12/22/2021
% Updated to enforce separability for control/estimator, SS/DJ - 01/13/2022
%


% % % 1. Check if a PIE has been specified
if ~exist('PIE','var')
    if exist('PIE_GUI','var')
        disp('No struct ''PIE'' specified, continuing with ''PIE_GUI''')
        PIE = PIE_GUI;
    elseif exist('PDE','var')
        disp('No struct ''PIE'' specified, continuing with ''PDE''')
        PIE = convert_PIETOOLS_PDE(PDE);
    elseif exist('PDE_GUI','var')
        disp('No struct ''PIE'' specified, continuing with ''PDE_GUI''')
        PIE = convert_PIETOOLS_PDE(PDE_GUI);
    else
        PIETOOLS_PDE_GUI
        error('Please specify a "PIE" struct (using the GUI), or call an example from the library, and run this script again')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 2. Determine what executive to run

% % Valid executives may be specified using one of the following
% stability=1;             % Test Lyapunov Stability
% stability_dual=1;        % Test Lyapunov Stability using an alternative Duality theorem
% Hinf_gain=1;             % Find a minimal bound on the Hing gain of the system
% Hinf_gain_dual=1;        % Find a minimal bound on the Hing gain of the system using an alternative duality theorem
% Hinf_control=1;          % Find a minimal bound on the Hinf optimal control problem
% Hinf_estimator=1;        % Find a minimal bound on the Hinf optimal observer problem

% % Check if an executive has been specified
exec = cell(0,0);
if exist('stability','var') && stability==1
    exec = [exec;'stability'];
end
if exist('stability_dual','var') && stability_dual==1
    exec = [exec;'stability_dual'];
end
if exist('Hinf_gain','var') && Hinf_gain==1
    exec = [exec;'Hinf_gain'];
end
if exist('Hinf_gain_dual','var') && Hinf_gain_dual==1
    exec = [exec;'Hinf_gain_dual'];
end
if exist('Hinf_estimator','var') && Hinf_estimator==1
    exec = [exec;'Hinf_estimator'];
end
if exist('Hinf_control','var') && Hinf_control==1
    exec = [exec;'Hinf_control'];
end

% % If no executive has been specified, let the user choose one
if isempty(exec)
    fprintf('\n What would you like to analyze/control? \n');
    msg = ['   Please input ''stability'', ''stability_dual'', ''Hinf_gain'', ''Hinf_gain_dual'', ''Hinf_estimator'', or ''Hinf_control'' \n ---> '];
    exec = input(msg,'s');
    exec = strrep(exec,'''','');
    exec = split(exec,[" ",","]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 3. Determine what settings to use

% % If no settings have been specified, allow the user to specify them
if ~exist('settings','var')
    fprintf('\n How hard should the solver work? \n');
    msg = ['   Please input ''extreme'', ''stripped'', ''light'', ''heavy'', ''veryheavy'', or ''custom'' \n ---> '];
    sttngs = input(msg,'s');
    sttngs = strrep(sttngs,'''','');    % Get rid of potential apostrophes
    settings = lpisettings(sttngs);     % Construct the settings
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 4. Run the executives

% % For each of the desired executives, define appropriate outputs, and run
for j=1:length(exec)
if strcmpi(exec{j},'stability')             % stability
    outval = '[prog_stability, P_stability]';
elseif strcmpi(exec{j},'stability_dual')    % stability_dual
    outval = '[prog_stability_d, P_stability_d]';
elseif strcmp(exec{j},'Hinf_gain')          % Hinf_gain
    outval = '[prog_Hinf_gain, P_Hinf_gain, Hinf_gain]';
elseif strcmp(exec{j},'Hinf_gain_dual')     % Hinf_gain_dual
    outval = '[prog_Hinf_gain_d, P_Hinf_gain_d, Hinf_gain_dual]';
elseif contains(exec{j},'control')          % Hinf_control
    outval = '[prog_control, K_control, Hinf_gain_control, P_control, Z_control]';
elseif contains(exec{j},'estimator')        % Hinf_estimator
    outval = '[prog_estimator, L_estimator, Hinf_gain_estimator, P_estimator, Z_estimator]';
end
infun = ['PIETOOLS_',exec{j},'(PIE,settings)'];
evalin('base',[outval,'=',infun,';']);  % Run the executive
end

% % % Clean up...
clear exec sttngs j msg infun outval