function display(PIE,name)
% Display the fields of an "npie_struct" object 'PIE', with name 'name'.
%
% INPUT
% - PIE:    npie_struct class object defining a nonlinear PIE.
% - name:   char specifying the name of the input 'PIE' in the workspace.
%
% OUTPUT
% Displays the different fields comprising the npie_struct in an ordered
% manner.
%
% NOTES
% This is a "disp" function, not a "display" function. It only shows which
% fields the input PIE pertains, ordered in a particular manner, without
% showing e.g. the structure of the associated system or the specific
% variables, domain or parameters. The sizes are the total number of
% variables in the function spaces PIE.space, which are specified by
% "initialize" based on the operators.

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
% DJ, 10/29/2026: Initial coding, based on the display function of
%                   pie_struct;

try
if nargin==1
    name = 'ans';
end
% Times character (x)
ttimes = '\xD7';

% Operators in each row map to the same space (x, z, or y), and operators
% in each column from the same space (x, w, or u).
opnames = {'T', 'Tw', 'Tu';
           'f', 'Bw', 'Bu';
           'gz','Dzw','Dzu';
           'gy','Dyw','Dyu'};
row_type = {'x';'x';'z';'y'};
col_type = {'x','w','u'};

% % % Display the PIE structure
% Start with the name
fprintf(['\n',name,' = \n'])
disp('  <a href="matlab:helpPopup npie_struct">npie_struct</a> with properties:')

% Then display the dimension, vars, and domain.
dim = PIE.dim;
fprintf(['\n     dim: ',num2str(dim),';\n',...
           '    vars: [',num2str(dim),ttimes,'2 polynomial];\n',...
           '     dom: [',num2str(dim),ttimes,'2 double];\n']);

% Then display the parameters.
for r=1:size(opnames,1)
    fprintf('\n');
    for c=1:size(opnames,2)
        op = PIE.(opnames{r,c});
        nr = sum(PIE.space.(row_type{r}).dims);
        nc = sum(PIE.space.(col_type{c}).dims);
        if isnumeric(op) && isempty(op)
            op_class = 'unspecified';
        else
            op_class = class(op);
        end
        fprintf(['  %3s: [%d',ttimes,'%d %s];   '],opnames{r,c},nr,nc,op_class);
    end
end
fprintf('\n');

% Display the operators defining the algebraic constraint, if any.
for field = {'K','Kw','Ku'}
    obj = PIE.(field{1});
    if size(obj,1)>0
        fprintf('\n%s: [%s %s];\n',field{1},mat2str(size(obj)),class(obj));
    end
end

if ~isempty(PIE.misc) && ~(isstruct(PIE.misc) && isempty(fieldnames(PIE.misc)))
    if isa(PIE.misc,'struct')
        fprintf('\n    misc\n');
        fnames = fieldnames(PIE.misc);
        for jj=1:numel(fnames)
            obj = PIE.misc.(fnames{jj});
            fprintf(['     .',fnames{jj},': [%d',ttimes,'%d %s];\n'],size(obj,1),size(obj,2),class(obj));
        end
    else
        fprintf(['\n    misc: [%d',ttimes,'%d %s];\n'],size(PIE.misc,1),size(PIE.misc,2),class(PIE.misc));
    end
end

fprintf('\n');

catch ME
    disp(PIE);
end

end
