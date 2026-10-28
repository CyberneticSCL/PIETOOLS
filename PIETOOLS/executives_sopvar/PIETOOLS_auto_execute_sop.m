%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_auto_execute_sop.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script that runs the container executives selected by flags in the
% workspace on a PIE, and stores their outputs in the base workspace.
%
% INPUT (variables in the workspace)
% - PIE:      'pie_struct' object; if absent, PIE_GUI ('pie_struct'), or
%             PDE or PDE_GUI ('pde_struct', converted by
%             convert_PIETOOLS_PDE); if none exists the GUI opens;
% - flags:    each a double, 1 to run:
%               stability       PIETOOLS_stability_sop, form 'Q';
%               stability_dual  PIETOOLS_stability_sop, form 'Qd';
%               Hinf_gain       PIETOOLS_Hinf_gain_sop, form 'Q';
%               Hinf_gain_dual  PIETOOLS_Hinf_gain_sop, form 'Qd';
%               Hinf_estimator, Hinf_control, H2_estimator, H2_control,
%               well_posedness  the executive of that name;
%               H2_norm         PIETOOLS_H2_norm_sop, form 'o';
%               H2_norm_dual    PIETOOLS_H2_norm_sop, form 'c';
%             with no flag set, the script asks which to run;
% - settings: (optional) 'struct' of LPI settings (lpisettings,
%             lpisettings_sop); if absent, the script asks for a level of
%             lpisettings ('light' when the answer is empty or unknown).
%
% OUTPUT (variables in the base workspace), per executive run:
%   stability:       prog_stability, P_stability, info_stability;
%   stability_dual:  prog_stability_d, P_stability_d, info_stability_d;
%   Hinf_gain:       prog_Hinf_gain, P_Hinf_gain, Hinf_gain, info_Hinf_gain;
%   Hinf_gain_dual:  prog_Hinf_gain_d, P_Hinf_gain_d, Hinf_gain_dual,
%                    info_Hinf_gain_d;
%   Hinf_control:    prog_control, K_control, Hinf_gain_control, P_control,
%                    Z_control, info_control;
%   Hinf_estimator:  prog_estimator, L_estimator, Hinf_gain_estimator,
%                    P_estimator, Z_estimator, info_estimator;
%   H2_norm:         prog_H2_norm, W_H2_norm, H2_norm, R_H2_norm,
%                    Q_H2_norm, info_H2_norm;
%   H2_norm_dual:    prog_H2_norm_d, W_H2_norm_d, H2_norm_dual,
%                    R_H2_norm_d, Q_H2_norm_d, info_H2_norm_d;
%   H2_control:      prog_control, K_control, H2_norm_control, P_control,
%                    Z_control, W_control, info_control;
%   H2_estimator:    prog_estimator, L_estimator, H2_norm_estimator,
%                    P_estimator, Z_estimator, W_estimator, info_estimator;
%   well_posedness:  prog_well_posedness, P_well_posedness,
%                    R_well_posedness, omega_well_posedness,
%                    info_well_posedness;
% the outputs of each executive as its header states them.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Initial coding MMP, 10/08/2026 (the stock script of DJ, 01/12/2025, with
% the _sop names and the info outputs)
% MMP, 10/28/2026: The executives are the consolidated ones of 10/28/2026:
%                stability, stability_dual, Hinf_gain_dual, H2_norm and
%                H2_norm_dual call PIETOOLS_stability_sop,
%                PIETOOLS_Hinf_gain_sop and PIETOOLS_H2_norm_sop with the
%                form of the executive they called before
%                (PIETOOLS_PIE2PDEstability_sop, _dual_sop,
%                PIETOOLS_Hinf_gain_dual_sop, PIETOOLS_H2_norm_o_sop,
%                PIETOOLS_H2_norm_c_sop, removed). The call string is set in
%                each branch; the former single line built it from the
%                name. Header to the header protocol.

% The PIE: from the workspace, the GUI, or a PDE
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

% The executives to run, by the flags in the workspace
exec = cell(0,0);
if exist('stability','var') && stability==1
    exec = [exec;'PIE2PDEstability'];
end
if exist('stability_dual','var') && stability_dual==1
    exec = [exec;'PIE2PDEstability_dual'];
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
if exist('H2_norm','var') && H2_norm==1
    exec = [exec;'H2_norm_o'];
end
if exist('H2_norm_dual','var') && H2_norm_dual==1
    exec = [exec;'H2_norm_c'];
end
if exist('H2_estimator','var') && H2_estimator==1
    exec = [exec;'H2_estimator'];
end
if exist('H2_control','var') && H2_control==1
    exec = [exec;'H2_control'];
end
if exist('well_posedness','var') && well_posedness==1
    exec = [exec;'well_posedness'];
end
if isempty(exec)
    fprintf('\n What would you like to analyze/control? \n');
    msg = ['   Please input ''stability'', ''stability_dual'', ''Hinf_gain'', ''Hinf_gain_dual'','...
            ' ''Hinf_estimator'', ''Hinf_control'', ''H2_norm'', ''H2_norm_dual'','...
            ' ''H2_estimator'', ''H2_control'', or ''well_posedness'' \n ---> '];
    exec = input(msg,'s');
    exec = strrep(exec,'''','');
    exec = split(exec,[" ",","]);
end

% The settings: from the workspace, or a level of lpisettings
if ~exist('settings','var')
    fprintf('\n How hard should the solver work? \n');
    msg = ['   Please input ''extreme'', ''stripped'', ''light'' (default), ''heavy'', ''veryheavy'', or ''custom''\n ---> '];
    sttngs = input(msg,'s');
    sttngs = strrep(sttngs,'''','');
    if isempty(sttngs)
        disp(' ')
        disp(" Proceeding with 'light' settings.")
        sttngs = 'light';
    elseif ~strcmp(sttngs,'extreme') && ~strcmp(sttngs,'stripped') && ~strcmp(sttngs,'light')...
            && ~strcmp(sttngs,'heavy') && ~strcmp(sttngs,'veryheavy') && ~strcmp(sttngs,'custom')
        fprintf(2,"\n Unknown settings specified; proceeding with 'light' settings.\n")
        sttngs = 'light';
    end
    settings = lpisettings(sttngs);
end

% Run each executive, assigning its outputs in the base workspace
for j=1:length(exec)
if strcmpi(exec{j},'PIE2PDEstability')
    outval = '[prog_stability, P_stability, info_stability]';
    infun = 'PIETOOLS_stability_sop(PIE,settings,''Q'')';
    msg_out = '';
elseif strcmpi(exec{j},'PIE2PDEstability_dual')
    outval = '[prog_stability_d, P_stability_d, info_stability_d]';
    infun = 'PIETOOLS_stability_sop(PIE,settings,''Qd'')';
    msg_out = '';
elseif strcmp(exec{j},'Hinf_gain')
    outval = '[prog_Hinf_gain, P_Hinf_gain, Hinf_gain, info_Hinf_gain]';
    infun = 'PIETOOLS_Hinf_gain_sop(PIE,settings,''Q'')';
    msg_out = 'The upper bound on the H-infty norm has been saved as "Hinf_gain" (certificate and witness in "info_Hinf_gain").';
elseif strcmp(exec{j},'Hinf_gain_dual')
    outval = '[prog_Hinf_gain_d, P_Hinf_gain_d, Hinf_gain_dual, info_Hinf_gain_d]';
    infun = 'PIETOOLS_Hinf_gain_sop(PIE,settings,''Qd'')';
    msg_out = 'The upper bound on the (dual) H-infty norm has been saved as "Hinf_gain_dual".';
elseif contains(exec{j},'Hinf_control')
    outval = '[prog_control, K_control, Hinf_gain_control, P_control, Z_control, info_control]';
    infun = 'PIETOOLS_Hinf_control_sop(PIE,settings)';
    msg_out = 'The (optimal) feedback control gain operator has been saved as "K_control".';
elseif contains(exec{j},'Hinf_estimator')
    outval = '[prog_estimator, L_estimator, Hinf_gain_estimator, P_estimator, Z_estimator, info_estimator]';
    infun = 'PIETOOLS_Hinf_estimator_sop(PIE,settings)';
    msg_out = 'The (optimal) estimator operator has been saved as "L_estimator".';
elseif contains(exec{j},'H2_norm_o')  ||  strcmp(exec{j},'H2_norm')
    outval = '[prog_H2_norm, W_H2_norm, H2_norm, R_H2_norm, Q_H2_norm, info_H2_norm]';
    infun = 'PIETOOLS_H2_norm_sop(PIE,settings,''o'')';
    msg_out = 'The upper bound on the H-2 norm has been saved as "H2_norm".';
elseif contains(exec{j},'H2_norm_c') || strcmp(exec{j},'H2_norm_dual')
    outval = '[prog_H2_norm_d, W_H2_norm_d, H2_norm_dual, R_H2_norm_d, Q_H2_norm_d, info_H2_norm_d]';
    infun = 'PIETOOLS_H2_norm_sop(PIE,settings,''c'')';
    msg_out = 'The upper bound on the (dual) H-2 norm has been saved as "H2_norm_dual".';
elseif contains(exec{j},'H2_control')
    outval = '[prog_control, K_control, H2_norm_control, P_control, Z_control, W_control, info_control]';
    infun = 'PIETOOLS_H2_control_sop(PIE,settings)';
    msg_out = 'The (optimal) feedback control gain operator has been saved as "K_control".';
elseif contains(exec{j},'H2_estimator')
    outval = '[prog_estimator, L_estimator, H2_norm_estimator, P_estimator, Z_estimator, W_estimator, info_estimator]';
    infun = 'PIETOOLS_H2_estimator_sop(PIE,settings)';
    msg_out = 'The (optimal) estimator operator has been saved as "L_estimator".';
elseif strcmpi(exec{j},'well_posedness')
    outval = '[prog_well_posedness, P_well_posedness, R_well_posedness, omega_well_posedness, info_well_posedness]';
    infun = 'PIETOOLS_well_posedness_sop(PIE,settings)';
    msg_out = '';
else
    error(['Unknown executive ''',exec{j},'''.'])
end
evalin('base',[outval,'=',infun,';']);
disp(msg_out);
end
clear exec sttngs j msg infun outval msg_out
