function display(P,name,width_max,call_name)
% DISPLAY(P,NAME,WIDTH_MAX,CALL_NAME) displays the PI operator P with
% decision variables in the command window as
%   (NAME x)(s) = int ... int K(s,s') x(s') ds' ... ds' + ...,
% with the operator written out term by term, each term showing the
% integrals taken, with the dummy variables displayed as gamma. Called
% automatically when displaying an 'sdopvar' object without semicolon.
%
% INPUTS
% - P:      'sdopvar' object representing a PI operator;
% - name:   (optional) char, the name to display for the operator. Defaults
%           to the name of the variable passed to display, or to 'Pop' if
%           there is none. The name 'ans' is replaced by 'Pop';
% - width_max:  (optional) maximal number of characters of a line with
%           multiple terms. Defaults to the width of the command window, or
%           to 100 if this cannot be determined;
% - call_name:  (optional) char, the expression with which the user can
%           call the operator, used in the call to extract a kernel.
%           Defaults to name, also if that is ans, and is set
%           by the display of a container that has P as its only block,
%           e.g. 'Pop.C{1,1}';
%
% OUTPUTS
% Display of the operator in the command window, with the terms constructed
% by 'get_display', see there for the integrals. The terms are placed on
% the same line, unless the line would have more characters than fit in the
% command window, in which case the whole next term is moved to the next
% line, so that no term is split by the command window, with the sign
% aligned under the equality sign.
%
% NOTES
% - A kernel with at most 4 nonzero coefficients is displayed explicitly,
%   if the operator is scalar-valued, with the decision variables displayed
%   by their names, as d1*s*g for the coefficient d1 of the monomial s*g.
%   For a matrix-valued operator, each element of the output is displayed
%   explicitly, (NAME x)_i(s) = ..., if each of its kernels has at most 4
%   nonzero coefficients per element, and the element has at most 100
%   characters; the elements of a vector-valued function x are named x_1,
%   x_2, .... Other kernels K_i are displayed by name. The kernel
%   K_i can be extracted by calling P.kernels{i}, a 'dpvar' object in the
%   variables, the dummy variables, and the decision variables P.Zd.
% - The function that the operator acts on is named x, unless a variable is
%   also named x (ignoring a number at the end of the name), in which case
%   it is named u, w, or v, the first that is not the name of a variable.
% - An identity operator is displayed as (NAME x)(s) = x(s), with the
%   spatial variables of the operator, and as NAME x = x if it acts on R^n.
% - The cost does not depend on the number of decision variables: only the
%   nonzeros of the parameters are accessed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - display
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
% DJ, 11/06/2026: Initial coding;
% DJ, 11/06/2026: Extract the kernels that are displayed by name with
%                   'kernels', instead of the coefficients;
% DJ, 11/06/2026: Get the subscripts from 'varname2display';

if nargin<2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'Pop';
    end
end
if nargin<4 || isempty(call_name)
    call_name = name;                   % the name that the user can call
end
if strcmp(name,'ans')
    % Display an unnamed object under a default name
    name = 'Pop';
    fprintf('\n  ans = %s, where\n',name);
end
if nargin<3
    width_max = [];
end
if isempty(P.vars) || prod(P.dims)==0
    fprintf('\n  %s is an empty operator, with input dimension %d and output dimension %d.\n\n',name,P.dims(2),P.dims(1));
    return
end

disp_struct = get_display(P);
if ~isempty(disp_struct.rows)
    out_names = cellfun(@varname2display,reshape(P.vars.out,1,[]),'UniformOutput',false);
    fprintf('\n');
    for i=1:numel(disp_struct.rows)
        lhs = [name,' ',disp_struct.xname];
        if ~isempty(out_names) || numel(disp_struct.rows)>1
            lhs = ['(',lhs,')'];
        end
        if numel(disp_struct.rows)>1
            [~,sub_i] = varname2display(sprintf('%d',i));
            lhs = [lhs,sub_i];
        end
        if ~isempty(out_names)
            lhs = [lhs,'(',strjoin(out_names,','),')'];
        end
        print_terms([lhs,' ='],disp_struct.rows(i).terms,disp_struct.rows(i).signs,width_max);
    end
    fprintf('\n');
    return
end
if isempty(P.vars.out)
    lhs = [name,' ',disp_struct.xname,' ='];
else
    out_names = cellfun(@varname2display,reshape(P.vars.out,1,[]),'UniformOutput',false);
    lhs = ['(',name,' ',disp_struct.xname,')(',strjoin(out_names,','),') ='];
end
fprintf('\n');
print_terms(lhs,disp_struct.terms,disp_struct.signs,width_max);

% Explain how to extract the kernels that are not displayed
if disp_struct.use_K
    if isscalar(P.params.A)
        fprintf('\n  Call "%s.kernels{1}" to extract kernel K.\n',call_name);
    else
        [~,sub_i] = varname2display('i');
        fprintf('\n  Call "%s.kernels{i}" to extract kernel K%s.\n',call_name,sub_i);
    end
end
fprintf('\n');

end
