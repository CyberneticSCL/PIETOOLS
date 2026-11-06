function print_terms(lhs,terms,signs,width_max)
% PRINT_TERMS(LHS,TERMS,SIGNS,WIDTH_MAX) prints the equation
%   LHS term_1 + term_2 - term_3 ...
% in the command window, with the terms on one line as long as they fit in
% the command window, and the whole next term on a new line otherwise.
%
% INPUTS
% - lhs:        char, the left-hand side of the equation, including the
%               equality sign, e.g. '(P x)(s) =';
% - terms:      1 x n cell of char, the terms of the sum, without their
%               sign;
% - signs:      (optional) 1 x n array of 1 and -1, the sign of each term.
%               Defaults to all 1;
% - width_max:  (optional) maximal number of characters of a line with
%               multiple terms. Defaults to the width of the command window
%               minus one, or to 100 if this cannot be determined;
%
% OUTPUTS
% The equation is printed. The terms are only moved to a new line as a
% whole, so that the command window does not split a term, with the sign
% aligned under the equality sign. A margin of one character is left, which
% is not used by the command window. Without terms, the right-hand side is
% 0.
%
% NOTES
% - The terms are used by the 'display' methods of the classes of
%   operators, which construct them with their 'get_display' methods.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - print_terms
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

% Default settings
width_default = 100;    % maximal width of a line, if the width of the command window
                        % cannot be determined
width_min = 40;         % minimal width of a line

if nargin<3 || isempty(signs)
    signs = ones(1,numel(terms));
end
if nargin<4 || isempty(width_max)
    width_max = width_default;
    try
        cmd_win_size = get(0,'CommandWindowSize');
        width_max = max(cmd_win_size(1)-1,width_min);
    catch
    end
end
indent_str = repmat(' ',1,2+numel(lhs)-nnz(lhs==char(770))-1);     % the combining circumflex takes no space
sign_chars = '+-';                  % sign_chars(1+(signs(j)<0)) is the sign of term j
if isempty(terms)
    fprintf('  %s 0\n',lhs);
    return
end
line_str = ['  ',lhs,' ',repmat('- ',1,signs(1)<0),terms{1}];
for j=2:numel(terms)
    term_sign = sign_chars(1+(signs(j)<0));
    if numel(line_str)+3+numel(terms{j})>width_max
        fprintf('%s\n',line_str);
        line_str = [indent_str,term_sign,' ',terms{j}];
    else
        line_str = [line_str,' ',term_sign,' ',terms{j}];
    end
end
fprintf('%s\n',line_str);

end
