function varList = getVarList(list, type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getPDEstructure.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes all variables list PDE, ODE, Inputs and Outputs and returns 
% the type specified by type.
% 
% INPUTS:
% list: GUI uitable structure with list of all variables
% type: 'distributed', 'finite', 'output', 'input', 'pde', 'ode',
% 'finiteRHS', 'allRHS', 'all' default is 'all'
% 
% OUTPUTS:
% varList: cell array of strings, list of values specified by type

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications

if ~isa(list, 'matlab.ui.control.Table')
    error("Input to getVarList must be a uitable");
end

% first get all variables out and the segregation
allVars = get(list,'Data');
headerLoc = get(list,'UserData');
headerLoc = headerLoc.headerLoc;

PDEvarlist = allVars(headerLoc(1,1)+1:headerLoc(2,1)-1,:);
ODEvarlist = allVars(headerLoc(2,1)+1:headerLoc(3,1)-1,:);
Disturbvarlist = allVars(headerLoc(3,1)+1:headerLoc(4,1)-1,:);
Inputvarlist = allVars(headerLoc(4,1)+1:headerLoc(5,1)-1,:);
RegOutvarlist = allVars(headerLoc(5,1)+1:headerLoc(6,1)-1,:);
ObsOutvarlist = allVars(headerLoc(6,1)+1:end,:);

if strcmp(type,'pde')|| strcmp(type,'distributed')
    varList = PDEvarlist;
elseif strcmp(type, 'ode')
    varList = ODEvarlist;
elseif strcmp(type, 'regulated')
    varList = RegOutvarlist;
elseif strcmp(type, 'observed')
    varList = ObsOutvarlist;
elseif strcmp(type, 'output')
    varList = [RegOutvarlist; ObsOutvarlist];
elseif strcmp(type, 'disturb')
    varList = Disturbvarlist;
elseif strcmp(type, 'control')
    varList = Inputvarlist;
elseif strcmp(type, 'input')
    varList = [Disturbvarlist; Inputvarlist];
elseif strcmp(type, 'finiteRHS')
    varList = [ODEvarlist; Disturbvarlist; Inputvarlist];
elseif strcmp(type, 'finite')
    varList = [ODEvarlist; Disturbvarlist; Inputvarlist; RegOutvarlist; ObsOutvarlist];
elseif strcmp(type, 'allRHS')
    varList = [PDEvarlist;ODEvarlist; Disturbvarlist; Inputvarlist];
else
    varList = [PDEvarlist;ODEvarlist; Disturbvarlist; Inputvarlist; RegOutvarlist; ObsOutvarlist];
end