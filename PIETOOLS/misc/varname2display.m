function [str,lo,hi] = varname2display(name,copy)
% [STR,LO,HI] = VARNAME2DISPLAY(NAME,COPY) returns the name NAME of a
% variable as it is displayed for objects of the operator classes (e.g.
% 'sopvar', 'intop', 'tensopvar', 'polyopvar') and for PIEs, and the forms
% of this name used as the limit of an integral, so that the display of
% names follows from this one function.
%
% INPUTS
% - name:   char, the name of a spatial variable, a state variable, a dummy
%           variable, or a number, e.g. 's1', 'x1_pie', 's_dum2', '0.5',
%           or an index, e.g. '2' or 'i';
% - copy:   (optional) positive integer, the number of the copy of the
%           variable, if it appears multiple times;
%
% OUTPUTS
% - str:    char, the name to display. A name with _pie (a state variable
%           of a PIE) gets a hat, displayed after the name, e.g. 'x_pie' as
%           x with a hat. Then follow the subscripts, in order,
%           - the number at the end of the name, e.g. 's1' is displayed as
%             s with the subscript 1, and 'x1_pie' as x with a hat and the
%             subscript 1;
%           - the number of the copy, separated by a comma, e.g. the copy 2
%             of the variable 's1' is displayed as s with the subscripts
%             1,2.
%           A name ending with _dum or _dum<k> is a dummy variable, which
%           is displayed as gamma, with the subscript k if given, e.g.
%           's_dum2' as gamma with the subscript 2. The rest of the name is
%           not used, so that the caller numbers the dummy variables. A
%           name starting with a digit, a sign or a point is a number,
%           and is displayed as is;
% - lo:     char, the name as a lower limit of an integral, which precedes
%           the integral symbol: the name as subscripts, if all of its
%           displayed characters have one (the digits, the minus, s, i, j
%           and gamma), e.g. for 's1' and 's_dum2', and {str}_ otherwise,
%           e.g. for 0.5 or for a copy. For an index, e.g. '2' or 'i', this
%           is the subscript to be used in the display of an operator,
%           e.g. K_2. The name may also be given as it is displayed, e.g.
%           the str of 's1', which gives the same lo;
% - hi:     char, the name as an upper limit of an integral, which follows
%           the integral symbol: as lo, but with superscripts, and ^{str}
%           otherwise. For an index, the superscript, e.g. R^2;
%
% NOTES
% - The hat is the combining circumflex U+0302, which is displayed over the
%   last character of the name, e.g. the x in x1_pie. It takes up a
%   character in the string, but no space in the display.
% - A name with _pie is a name in which _pie is preceded by at least one
%   character, and followed only by an (optional) number, e.g. x_pie or
%   x_pie2, not x_piece.
% - If a name displayed with a subscript is given, e.g. the str of x1, the
%   subscript of a copy is separated from it by a comma.
% - Characters are specified by their code, to keep this file ASCII.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - varname2display
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
% DJ, 11/06/2026: Add dummy variables and numbers, and the forms for the
%                   limits of an integral;

sub_dig = char(8320 + (0:9));
sbs = @(n) sub_dig(n-'0'+1);

% Name, hat and index
pie = '';       idx = '';
if ~isempty(regexp(name,'^[-+.0-9]','once'))
    base = name;                                % number
else
    tok = regexp(name,'^.+_dum(\d*)$','tokens','once');
    if ~isempty(tok)
        base = char(947);                       % dummy variable
        idx = sbs(tok{1});
    else
        base = name;
        tok = regexp(base,'^(.+?)_pie_?(\d*)$','tokens','once');
        if ~isempty(tok)
            base = [tok{1},tok{2}];
            pie = char(770);
        end
        tok = regexp(base,'^(.*[^0-9_])_?(\d+)$','tokens','once');
        if ~isempty(tok)
            base = tok{1};
            idx = sbs(tok{2});
        end
    end
end
cp = '';
if nargin>=2 && ~isempty(copy)
    cp = sbs(sprintf('%d',copy));
end

% Combine the name, the hat, and the subscripts separated by commas, with a
% comma after a name that already has a subscript
sub_strs = {idx,cp};
sub_strs = sub_strs(~cellfun(@isempty,sub_strs));
sep = repmat(',',1,~isempty(sub_strs) && ~isempty(base) && ismember(base(end),sub_dig));
str = [base,pie,sep,strjoin(sub_strs,',')];

% The name as a limit of an integral: as subscript or superscript if all
% characters of the name have one (also a subscript digit), and in braces
% otherwise, as for the display of a PDE
glyphs = ['0123456789',sub_dig,'-s',char([947,105,106])];
subs = char([8320:8329,8320:8329,8331,8347,7527,7522,11388]);
sups = char([8304,185,178,179,8308:8313,8304,185,178,179,8308:8313,8315,738,7518,8305,690]);
[is_glyph,pos] = ismember(str,glyphs);
if all(is_glyph)
    lo = subs(pos);
    hi = sups(pos);
else
    lo = ['{',str,'}_'];
    hi = ['^{',str,'}'];
end

end
