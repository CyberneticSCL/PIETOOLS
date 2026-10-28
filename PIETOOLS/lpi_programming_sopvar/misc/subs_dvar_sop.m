function Xsol = subs_dvar_sop(X,dnames,dvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XSOL = SUBS_DVAR_SOP(X,DNAMES,DVALS) substitutes the value DVALS(k) for
% the decision variable named DNAMES{k} in X, for every k, without a solved
% program: XSOL is the solved value of X in a program whose decision
% variables are DNAMES, with the values DVALS.
%
% INPUT
% - X:      the object to substitute in: an 'sdopvar', 'cdopvar', 'dopvar',
%           'dopvar2d' or 'dpvar' object; an 'opvar' or 'opvar2d' object,
%           whose fields may hold 'dpvar' objects; a 'sopvar', 'copvar',
%           'polynomial' or 'double' object; a 'char' decision variable
%           name; or a 'cell' of these (see NOTES);
% - dnames: 'cellstr' or string array of distinct decision variable names,
%           in any order; names X does not hold are not used;
% - dvals:  numeric array with numel(dnames) entries, dvals(k) the value of
%           dnames{k}.
%
% OUTPUT
% - Xsol:   X with every decision variable named in DNAMES replaced by its
%           value; the decision variables of X not named stay decision
%           variables. An 'sdopvar' or 'cdopvar' X gives a 'sopvar' or
%           'copvar' object when DNAMES names all of its decision
%           variables, else an 'sdopvar' or 'cdopvar' object on the
%           others; a 'dopvar' or 'dopvar2d' X with all of them named gives
%           an 'opvar' or 'opvar2d' object; a 'dpvar' X gives a
%           'polynomial' object (all named) or a 'dpvar' object on the
%           others; the 'dpvar' fields of an 'opvar' or 'opvar2d' X are
%           substituted; a 'char' name gives a 'polynomial' object, its
%           value when it is named; a 'cell' is substituted entry by
%           entry; a 'sopvar', 'copvar', 'polynomial' or 'double' X, and an
%           'opvar' or 'opvar2d' X without 'dpvar' fields, are returned as
%           they are.
%
% NOTES
% A 'cell' X that holds no 'sopvar', 'sdopvar', 'copvar' or 'cdopvar'
% object, at any depth, may hold only 'dpvar', 'polynomial', 'double' and
% 'char' objects and cells of them, and its 'polynomial' entries are
% substituted. DNAMES not a 'cellstr' or string array, or DVALS not numeric
% or with a number of entries other than numel(DNAMES), is an error.
%
% See also LPIGETSOL_SOP, GETSOL_LPIVAR_SOP, SOSGETSOL.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - subs_dvar_sop
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
% Initial coding MMP, 09/29/2026. Decision-vector substitution without a
%                solved program, for certification of candidate points;
%                the low-rank certifiers (pielr_evalop/pielr_certify) had
%                it for dopvar/dopvar2d/dpvar only.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The deleted header text was the purpose (the
%                certification of candidate points), the rationale of
%                running lpigetsol_sop, and the cost by reference to
%                getsol_lpivar_sop and lpigetsol.

% DNAMES as a 'cellstr'; DVALS numeric, one entry per name
if isstring(dnames),    dnames = cellstr(dnames);   end
if ~iscellstr(dnames) %#ok<ISCLSTR>
    error('subs_dvar_sop:badNames','DNAMES must be a cellstr or string array of decision variable names.')
end
if ~isnumeric(dvals) || numel(dvals)~=numel(dnames)
    error('subs_dvar_sop:badVals','DVALS must be numeric with one entry per name (%d names, %d values).',...
        numel(dnames),numel(dvals))
end
% The fields lpigetsol_sop, getsol_lpivar_sop, sosgetsol and getsol_lpivar
% read, and no other; 'info' only marks the structure as solved.
prog = struct('decvartable',{dnames(:)},...
              'solinfo',struct('RRx',double(dvals(:)),'info',struct('source','subs_dvar_sop')));
% The solved value of X for that program
Xsol = lpigetsol_sop(prog,X);

end
