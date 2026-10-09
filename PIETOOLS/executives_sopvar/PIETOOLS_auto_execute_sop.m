%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS_auto_execute_sop.m     PIETOOLS 2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The container-path counterpart of PIETOOLS_auto_execute: a script that
% reads the same flags from the workspace (stability, stability_dual,
% Hinf_gain, Hinf_gain_dual, Hinf_estimator, Hinf_control, H2_norm,
% H2_norm_dual, H2_estimator, H2_control, well_posedness), the struct 'PIE'
% (or PIE_GUI, PDE, PDE_GUI) and 'settings' (or asks for a preset), and
% calls the PIETOOLS_*_sop executives of this folder with the stock output
% names plus one 'info_*' struct per executive (certificate status, degree
% loop history, numerical witness, dual kernel, SDP shape).
%
% Initial coding MMP, 10/08/2026 (the stock script of DJ, 01/12/2025, with
% the _sop names and the info outputs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% The settings
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

% Run each executive, assigning the stock output names and an info struct
for j=1:length(exec)
if strcmpi(exec{j},'PIE2PDEstability')
    outval = '[prog_stability, P_stability, info_stability]';
    msg_out = '';
elseif strcmpi(exec{j},'PIE2PDEstability_dual')
    outval = '[prog_stability_d, P_stability_d, info_stability_d]';
    msg_out = '';
elseif strcmp(exec{j},'Hinf_gain')
    outval = '[prog_Hinf_gain, P_Hinf_gain, Hinf_gain, info_Hinf_gain]';
    msg_out = 'The upper bound on the H-infty norm has been saved as "Hinf_gain" (certificate and witness in "info_Hinf_gain").';
elseif strcmp(exec{j},'Hinf_gain_dual')
    outval = '[prog_Hinf_gain_d, P_Hinf_gain_d, Hinf_gain_dual, info_Hinf_gain_d]';
    msg_out = 'The upper bound on the (dual) H-infty norm has been saved as "Hinf_gain_dual".';
elseif contains(exec{j},'Hinf_control')
    outval = '[prog_control, K_control, Hinf_gain_control, P_control, Z_control, info_control]';
    msg_out = 'The (optimal) feedback control gain operator has been saved as "K_control".';
elseif contains(exec{j},'Hinf_estimator')
    outval = '[prog_estimator, L_estimator, Hinf_gain_estimator, P_estimator, Z_estimator, info_estimator]';
    msg_out = 'The (optimal) estimator operator has been saved as "L_estimator".';
elseif contains(exec{j},'H2_norm_o')  ||  strcmp(exec{j},'H2_norm')
    exec{j} = 'H2_norm_o';
    outval = '[prog_H2_norm, W_H2_norm, H2_norm, R_H2_norm, Q_H2_norm, info_H2_norm]';
    msg_out = 'The upper bound on the H-2 norm has been saved as "H2_norm".';
elseif contains(exec{j},'H2_norm_c') || strcmp(exec{j},'H2_norm_dual')
    exec{j} = 'H2_norm_c';
    outval = '[prog_H2_norm_d, W_H2_norm_d, H2_norm_dual, R_H2_norm_d, Q_H2_norm_d, info_H2_norm_d]';
    msg_out = 'The upper bound on the (dual) H-2 norm has been saved as "H2_norm_dual".';
elseif contains(exec{j},'H2_control')
    outval = '[prog_control, K_control, H2_norm_control, P_control, Z_control, W_control, info_control]';
    msg_out = 'The (optimal) feedback control gain operator has been saved as "K_control".';
elseif contains(exec{j},'H2_estimator')
    outval = '[prog_estimator, L_estimator, H2_norm_estimator, P_estimator, Z_estimator, W_estimator, info_estimator]';
    msg_out = 'The (optimal) estimator operator has been saved as "L_estimator".';
elseif strcmpi(exec{j},'well_posedness')
    outval = '[prog_well_posedness, P_well_posedness, R_well_posedness, omega_well_posedness, info_well_posedness]';
    msg_out = '';
else
    error(['Unknown executive ''',exec{j},'''.'])
end
infun = ['PIETOOLS_',exec{j},'_sop(PIE,settings)'];
evalin('base',[outval,'=',infun,';']);
disp(msg_out);
end
clear exec sttngs j msg infun outval msg_out
