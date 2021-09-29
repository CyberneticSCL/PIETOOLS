%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIESIM_initial_setup.m    PIETOOLS 2021d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds time-derivatives of the user-derfined inputs and disturbances and
% adds them to the field 'uinput'
%
% Input:
% uinput - user-defined boundary inputs and disturbances
%
% Output:
% uinput - temporal derivatives of the user-defined boundary inputs and disturbances are added to the 
% field 'input' - these are needed for temporal integration of PIE equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding YP  - 9_18_2021

function uinput=PIESIM_initial_setup(uinput)
syms st;


 if isfield(uinput,'w')
             uinput.wdot=diff(uinput.w,st);
 end
 if isfield(uinput,'u')
             uinput.udot=diff(uinput.u,st);
 end
