function check_opts(opts,fname,defaults)
% CHECK_OPTS(OPTS,FNAME,DEFAULTS) warns about fields of the options OPTS of
% the function FNAME that are not recognized (e.g. misspelled fields),
% which would otherwise be ignored without notice.
%
% INPUTS
% - opts:       'struct' of options; anything else is not checked;
% - fname:      'char', the name of the function, used in the warning;
% - defaults:   (optional) 'struct' with all recognized fields, whose
%               'struct' fields list the recognized subfields, e.g. the
%               output of PIESOS_settings. Defaults to the options of the
%               PIESOS programming functions (piesos_ineq, piesos_sosvar,
%               piesos_quadvar, piesos_LF, ...), also within opts.L and
%               opts.R;
%
% NOTES
% - Options are passed on from one PIESOS function to another: only the
%   outermost function called by the user warns.
% - The warning has identifier 'PIESOS:unknownOption'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - check_opts
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

if ~isstruct(opts) || isempty(fieldnames(opts))
    return
end
% Only the outermost PIESOS function warns (frames above its caller)
st = dbstack(2);
if any(strncmpi({st.name},'piesos_',7) | strncmp({st.name},'monopmat',8))
    return
end

% Recognized fields, and the fields with recognized subfields
if nargin<3
    known = {'deg','degs','derivative','exclude','intmat','L','otype','pos','psatz','R','sep','sym','tol','unique'};
    subs = {'L','R'};
else
    known = fieldnames(defaults);
    subs = known(structfun(@isstruct,defaults));
end
bad = setdiff(fieldnames(opts),known);
for k=1:numel(subs)
    if isfield(opts,subs{k}) && isstruct(opts.(subs{k}))
        if nargin<3
            known_k = known;
        else
            known_k = fieldnames(defaults.(subs{k}));
        end
        bad = [bad; strcat(subs{k},'.',setdiff(fieldnames(opts.(subs{k})),known_k))];     %#ok<AGROW>
    end
end
if ~isempty(bad)
    warning('PIESOS:unknownOption',"%s: unrecognized option field(s) %s, which are ignored.",fname,strjoin(bad,', '))
end

end
