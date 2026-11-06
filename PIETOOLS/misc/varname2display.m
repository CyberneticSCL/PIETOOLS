function str = varname2display(name,copy)
% STR = VARNAME2DISPLAY(NAME,COPY) returns the name NAME of a variable as
% it is displayed for objects of the operator classes (e.g. 'sopvar',
% 'intop', 'tensopvar', 'polyopvar') and for PIEs, with subscripts for the
% index in the name, the type of variable (a PIE state or a dummy
% variable), and the copy of the variable.
%
% INPUTS
% - name:   char, the name of a spatial variable, a state variable, or a
%           dummy variable, e.g. 's1', 'x1_pie', 's1_dum2';
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
%           - the number k for a name ending with _dum<k> (the kth dummy
%             variable of a kernel), separated from the previous subscripts
%             by a comma, and followed by a prime after all of the
%             subscripts, e.g. 's1_dum2' is displayed as s with the
%             subscripts 1,2 and the prime;
%           - the number of the copy, separated by a comma, e.g. the copy 2
%             of the variable 's1' is displayed as s with the subscripts
%             1,2.
%
% NOTES
% - The hat is the combining circumflex U+0302, which is displayed over the
%   last character of the name, e.g. the x in x1_pie. It takes up a
%   character in the string, but no space in the display.
% - A name with _pie is a name in which _pie is preceded by at least one
%   character, and followed only by an (optional) number, e.g. x_pie or
%   x_pie2, not x_piece.
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

sub_dig = char(8320 + (0:9));
pie_hat = char(770);
prime = char(8242);
sbs = @(n) sub_dig(n-'0'+1);

base = name;
dum = '';       is_dum = false;
pie = '';
tok = regexp(base,'^(.+)_dum(\d+)$','tokens','once');
if ~isempty(tok)
    base = tok{1};
    dum = sbs(tok{2});
    is_dum = true;
end
tok = regexp(base,'^(.+?)_pie_?(\d*)$','tokens','once');
if ~isempty(tok)
    base = [tok{1},tok{2}];
    pie = pie_hat;
end
idx = '';
tok = regexp(base,'^(.*[^0-9_])_?(\d+)$','tokens','once');
if ~isempty(tok)
    base = tok{1};
    idx = sbs(tok{2});
end
cp = '';
if nargin>=2 && ~isempty(copy)
    cp = sbs(sprintf('%d',copy));
end
% Combine the name, the hat, the subscripts separated by commas, and the prime
sub_strs = {idx,dum,cp};
sub_strs = sub_strs(~cellfun(@isempty,sub_strs));
str = [base,pie,strjoin(sub_strs,','),repmat(prime,1,is_dum)];

end
