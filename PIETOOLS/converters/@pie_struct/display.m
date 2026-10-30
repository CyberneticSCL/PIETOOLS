function display(PIE,name)
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
% PIETOOLS - display
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 06/25/2025: Add display of 'misc' field;
% DB, 10/27/2026: Display optional PI operators F/K in case the null
% space of the derivative operator is non-trivial.
% DJ, 10/30/2026: Update to new pie_struct structure (display the function
%                   spaces, and constraint operators if specified);

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
nx = sum(PIE.space.x.dims);                                                 % DJ, 10/30/2026
nw = sum(PIE.space.w.dims);
nu = sum(PIE.space.u.dims);
nz = sum(PIE.space.z.dims);
ny = sum(PIE.space.y.dims);

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
       
% DJ, 10/29/2026
% Operators are named by their input and output (Bw, Cz, ...), and are
% typically copvar objects.
if isa(PIE.T,'opvar')
    op = 'opvar';
elseif isa(PIE.T,'opvar2d')
    op = 'opvar2d';
else
    op = 'copvar';
end

% Then display the parameters.
fprintf(['\n       T: [',nx,ttimes,nx,' ',op,'];     ',...
                 'Tw: [',nx,ttimes,nw,' ',op,'];     ',...
                 'Tu: [',nx,ttimes,nu,' ',op,']; ']);

fprintf(['\n       A: [',nx,ttimes,nx,' ',op,'];     ',...
                 'Bw: [',nx,ttimes,nw,' ',op,'];     ',...
                 'Bu: [',nx,ttimes,nu,' ',op,']; ']);

fprintf(['\n      Cz: [',nz,ttimes,nx,' ',op,'];    ',...
                'Dzw: [',nz,ttimes,nw,' ',op,'];    ',...
                'Dzu: [',nz,ttimes,nu,' ',op,']; ']);

fprintf(['\n      Cy: [',ny,ttimes,nx,' ',op,'];    ',...
                'Dyw: [',ny,ttimes,nw,' ',op,'];    ',...
                'Dyu: [',ny,ttimes,nu,' ',op,']; ']);

fprintf('\n');

% Display the operators defining the algebraic constraint, if any.
nK = size(PIE.K,1);                                                         % DJ, 10/30/2026
if nK>0
    nK = num2str(nK);
    fprintf(['\n       K: [',nK,ttimes,nx,' ',op,'];     ',...
                     'Kw: [',nK,ttimes,nw,' ',op,'];     ',...
                     'Ku: [',nK,ttimes,nu,' ',op,']; \n']);
end

if ~isempty(PIE.misc) && ~(isstruct(PIE.misc) && isempty(fieldnames(PIE.misc)))
    if isa(PIE.misc,'struct')
        fprintf('\n    misc\n');
        fnames = fieldnames(PIE.misc);
        for jj=1:numel(fnames)
            fprintf(['     .',fnames{jj},': ']);
            obj = PIE.misc.(fnames{jj});
            obj_class = class(obj);
            nr = cell2mat(num_list(str2num(num2str(size(obj,1))')+1)');
            nc = cell2mat(num_list(str2num(num2str(size(obj,2))')+1)');
            fprintf(['[',nr,ttimes,nc,' ',obj_class,'];\n'])
        end
    else
        fprintf('\n    misc: ')
        obj = PIE.misc;
        obj_class = class(obj);
        nr = cell2mat(num_list(str2num(num2str(size(obj,1))')+1)');
        nc = cell2mat(num_list(str2num(num2str(size(obj,2))')+1)');
        fprintf(['[',nr,ttimes,nc,' ',obj_class,'];\n'])
    end
end

fprintf('\n');

catch ME
    disp(PIE);
end
