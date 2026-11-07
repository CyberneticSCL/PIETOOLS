function out = use_unicode(val)
% OUT = USE_UNICODE(VAL) sets, or returns, the setting that determines
% whether the displays of operators and PIEs use special symbols that are
% missing from many fonts, see 'varname2display'.
%
% INPUTS
% - val:    (optional) logical, true to use the special symbols: subscript
%           and superscript letters, the hat on the states and inputs of a
%           PIE, and the circled times for the tensor product. False to
%           display these with plain characters. The value is stored as a
%           preference, so that it also applies in later MATLAB sessions.
%           Without input, the current setting is returned;
%
% OUTPUTS
% - out:    logical, the current setting, false by default. Only returned
%           if there is no input, or if an output is requested;
%
% NOTES
% - If false, the hat on the states and inputs of a PIE is displayed as an
%   f, e.g. x_pie as xf, the limits of an integral as {s_1}_ and ^{s_1},
%   instead of as subscript and superscript letters, the derivative in
%   time as d_t, the tensor product as o, and a tensor power as ^{o2}. The
%   subscript digits, gamma, and the integral sign are used in both cases.
% - The setting is read once per MATLAB session. Call 'clear use_unicode'
%   to read the stored preference again.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - use_unicode
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

% Default setting
default_unicode = false;      % display with plain characters

persistent is_unicode
if isempty(is_unicode)
    is_unicode = default_unicode;
    try
        is_unicode = logical(getpref('PIETOOLS','use_unicode',default_unicode));
    catch
    end
end
if nargin>=1
    if ischar(val) && strcmp(val,'true')        % allow "use_unicode true"
        val = true;
    elseif ischar(val) && strcmp(val,'false')
        val = false;
    end
    if ~(islogical(val) || isnumeric(val)) || ~isscalar(val)
        error('The setting should be a logical scalar.')
    end
    is_unicode = logical(val);
    try
        setpref('PIETOOLS','use_unicode',is_unicode);
    catch
        warning('The setting could not be stored, and only holds for this session.')
    end
end
if nargout>0 || nargin==0
    out = is_unicode;
end

end
