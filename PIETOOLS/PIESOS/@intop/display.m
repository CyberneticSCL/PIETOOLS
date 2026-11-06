function display(Kop,name,width_max,call_name)
% DISPLAY(KOP,NAME,WIDTH_MAX,CALL_NAME) displays the functional KOP in the
% command window as
%   NAME x = int ... int K(s) x(s) ds ... ds + ...,
% with the integral operator written out term by term, each term showing
% the integrals taken, in the order in which they are taken. Called
% automatically when displaying an 'intop' object without semicolon.
%
% INPUTS
% - Kop:    'intop' object representing a functional;
% - name:   (optional) char, the name to display for the functional.
%           Defaults to the name of the variable passed to display, or to
%           'Kop' if there is none. The name 'ans' is replaced by 'Kop';
% - width_max:  (optional) maximal number of characters of a line with
%           multiple terms. Defaults to the width of the command window, or
%           to 100 if this cannot be determined;
% - call_name:  (optional) char, the expression with which the user can
%           call the functional, used in the call to extract a kernel.
%           Defaults to name, also if that is ans, and is set by the
%           display of a container that has Kop as its only block, e.g.
%           'Rop.ops{1,1}';
%
% OUTPUTS
% Display of the functional in the command window, with the terms
% constructed by 'get_display', see there for the integrals. The terms are
% placed on the same line, unless the line would have more characters than
% fit in the command window, in which case the whole next term is moved to
% the next line, so that no term is split by the command window, with the
% sign aligned under the equality sign.
%
% NOTES
% - The integrals involve only dummy variables, displayed as gamma, with
%   the number of the integrated copy as subscript if there are multiple,
%   see 'display_term'. The domains of the variables are those of Kop.dom.
% - A kernel with at most 4 nonzero coefficients is displayed explicitly,
%   if the functional is scalar-valued. Other kernels K_i are displayed by
%   name; the kernel K_i can be extracted by calling Kop.kernels{i}, which
%   returns the kernel defined by Kop.params{i}, in the output variables
%   and in dummy variables named as the variables in Kop.vars with _dum
%   appended.
% - The function that the functional acts on is named x.
% - The terms are numbered like the elements of Kop.params, so a term
%   with a zero kernel is not displayed, but does take an index.
% - The cost does not depend on the number of decision variables: only the
%   nonzeros of the parameters of kernels that are displayed explicitly
%   are accessed.
% - An identity operator is displayed as NAME x = x.

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
% DJ, 11/06/2026: Get the subscripts from 'varname2display';

if nargin<2 || isempty(name)
    name = inputname(1);
    if isempty(name)
        name = 'Kop';
    end
end
if nargin<4 || isempty(call_name)
    call_name = name;                   % the name that the user can call
end
if strcmp(name,'ans')
    % Display an unnamed object under a default name
    name = 'Kop';
    fprintf('\n  ans = %s, where\n',name);
end

% An improperly declared functional cannot be displayed as an operator
if ~Kop.chkval && ~isvalid(Kop)
    builtin('disp',Kop);
    return
end
[m,n] = size(Kop);
if m*n==0
    fprintf('\n  %s is an empty functional, with input dimension %d and output dimension %d.\n\n',name,n,m);
    return
end

if nargin<3
    width_max = [];
end
disp_struct = get_display(Kop);
lhs = [name,' ',disp_struct.xname,' ='];
fprintf('\n');
print_terms(lhs,disp_struct.terms,disp_struct.signs,width_max);

if disp_struct.use_K
    if isscalar(Kop.params)
        fprintf('\n  Call "%s.kernels{1}" to extract kernel K.\n',call_name);
    else
        [~,sub_i] = varname2display('i');
        fprintf('\n  Call "%s.kernels{i}" to extract kernel K%s.\n',call_name,sub_i);
    end
end
fprintf('\n');

end
