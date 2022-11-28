function display(PIE,name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the fields of a "pie_struct" object 'PIE', with name 'name'.
%
% INPUT
% - PIE:    pie_struct class object defining a PIE.
% - name:   char specifying the name of the input 'PIE' in the workspace.
%
% OUTPUT
% Displays the different fields comprising the pie_struct in an ordered
% manner.
%
% NOTES
% This is a "disp" function, not a "display" function. It only shows which
% fields the input PIE pertains, ordered in a particular manner, without
% showing e.g. the structure of the associated system or the specific
% variables, domain or parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - disp_pie
%
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
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding DJ - 08/30/2022

% % % Set up a list of characters that we need for the display
% Name of the PIE structure.

try 
if nargin==1
    name = 'ans';
end
% Times character (x)
ttimes = '\xD7';
% Numbers 0 through 9
num_list = {'\x30';
           '\x31';
           '\x32';
           '\x33';
           '\x34';
           '\x35';
           '\x36';
           '\x37';
           '\x38';
           '\x39'};


% % % Extract information regarding the size of the problem.
% Check the spatial dimension of the problem.
dim = PIE.dim;

% Check the sizes of the states, inputs and outputs.
nx = [size(PIE.T,1); size(PIE.Tw,1); size(PIE.Tu,1); size(PIE.T,2);
      size(PIE.A,1); size(PIE.B1,1); size(PIE.B2,1); size(PIE.A,2);
      size(PIE.C1,2); size(PIE.C2,2)];
nw = [size(PIE.Tw,2); size(PIE.B1,2); size(PIE.D11,2); size(PIE.D21,2)];
nu = [size(PIE.Tu,2); size(PIE.B2,2); size(PIE.D12,2); size(PIE.D22,2)];
nz = [size(PIE.C1,1); size(PIE.D11,1); size(PIE.D12,1)];
ny = [size(PIE.C2,1); size(PIE.D21,1); size(PIE.D22,1)];

if any(nx(2:end)-nx(1:end-1))
    error('The PIE is not correctly specified: The dimensions of the fundamental state as suggested by the operators do not match.')
else
    nx = nx(1);
end
if any(nw(2:end)-nw(1:end-1))
    error('The PIE is not correctly specified: The dimensions of the exogenous input as suggested by the operators do not match.')
else
    nw = nw(1);
end
if any(nu(2:end)-nu(1:end-1))
    error('The PIE is not correctly specified: The dimensions of the actuator input as suggested by the operators do not match.')
else
    nu = nu(1);
end
if any(nz(2:end)-nz(1:end-1))
    error('The PIE is not correctly specified: The dimensions of the regulated output as suggested by the operators do not match.')
else
    nz = nz(1);
end
if any(ny(2:end)-ny(1:end-1))
    error('The PIE is not correctly specified: The dimensions of the observed output as suggested by the operators do not match.')
else
    ny = ny(1);
end

% Convert digits to char objects
nx = cell2mat(num_list(str2num(num2str(nx)')+1)');
nw = cell2mat(num_list(str2num(num2str(nw)')+1)');
nu = cell2mat(num_list(str2num(num2str(nu)')+1)');
nz = cell2mat(num_list(str2num(num2str(nz)')+1)');
ny = cell2mat(num_list(str2num(num2str(ny)')+1)');


% % % Display the PIE structure
% Start with the name
fprintf(['\n',name,' = \n'])
%disp('<a href="matlab: opentoline(which(''pie_struct.m''),1)">pie_struct</a> with properties:')
disp('  <a href="matlab:helpPopup pie_struct">pie_struct</a> with properties:')        

% Then display the dimension, vars, and domain.
fprintf(['\n     dim: ',num2str(dim),';\n',...
           '    vars: [',num2str(dim),ttimes,num_list{2+1},' polynomial];\n',...
           '     dom: [',num2str(dim),ttimes,num_list{2+1},' double];\n']);
       
if isa(PIE.T,'opvar')
    op = 'opvar';
elseif isa(PIE.T,'opvar2d')
    op = 'opvar2d';
end

% Then display the parameters.
fprintf(['\n       T: [',nx,ttimes,nx,' ',op,'];     ',...
                 'Tw: [',nx,ttimes,nw,' ',op,'];     ',...
                 'Tu: [',nx,ttimes,nw,' ',op,']; ']);

fprintf(['\n       A: [',nx,ttimes,nx,' ',op,'];     ',...
                 'B1: [',nx,ttimes,nw,' ',op,'];     ',...
                 'B2: [',nx,ttimes,nu,' ',op,']; ']);

fprintf(['\n      C1: [',nz,ttimes,nx,' ',op,'];    ',...
                'D11: [',nz,ttimes,nw,' ',op,'];    ',...
                'D12: [',nz,ttimes,nu,' ',op,']; ']);

fprintf(['\n      C2: [',ny,ttimes,nx,' ',op,'];    ',...
                'D21: [',ny,ttimes,nw,' ',op,'];    ',...
                'D22: [',ny,ttimes,nu,' ',op,']; ']);

fprintf('\n');

catch ME
    disp(PIE);
end
