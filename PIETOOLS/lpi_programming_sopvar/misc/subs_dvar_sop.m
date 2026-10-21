function Xsol = subs_dvar_sop(X,dnames,dvals)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% XSOL = SUBS_DVAR_SOP(X,DNAMES,DVALS) substitutes given values for named
% decision variables in an LPI object X of either operator family, with no
% solved program: decision variable DNAMES{k} takes the value DVALS(k).
%
% This is the entry point for certification, where the point is not a
% solver's answer: a candidate from a low-rank or first-order method, a
% rounded point, or a perturbed one. It is 'lpigetsol_sop' applied to a
% program structure holding only DNAMES as decvartable and DVALS as
% solinfo.RRx, so both entry points run the one implementation, and every
% class 'lpigetsol_sop' accepts is accepted here with the same result.
%
% INPUT
% - X:      object of any class accepted by 'lpigetsol_sop': 'sdopvar',
%           'cdopvar', 'dopvar', 'dopvar2d', 'dpvar', a fixed operator or
%           matrix, or a cell of these;
% - dnames: cellstr (or string array) of decision variable names, in any
%           order, UNIQUE. Names not in X are ignored;
% - dvals:  numeric vector, dvals(k) the value of dnames{k}.
%
% OUTPUT
% - Xsol:   X with every listed decision variable replaced by its value.
%           Variables of X not in DNAMES stay decision variables, so a
%           subset of the variables can be fixed (e.g. only gamma).
%
% COST: as 'getsol_lpivar_sop' with N = numel(dnames); for the legacy
% classes, as 'lpigetsol'.
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
Xsol = lpigetsol_sop(prog,X);

end
