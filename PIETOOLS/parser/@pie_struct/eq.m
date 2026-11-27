function logval=eq(objA, objB)
% LOGVAL = EQ(OBJA,OBJB) or OBJA==OBJB tests whether two "pie_struct"
% objects represent the same PIE.
%
% INPUT
% - objA, objB:    pie_struct class objects defining PIEs;
%
% OUTPUT
% - logval: true if the PIEs involve the same spatial variables and domain,
%           and each of the operators is the same, false otherwise.
%
% NOTES
% Operators of different types, or acting on different function spaces,
% are not equal.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - eq
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% Initial coding PIETOOLS Team (PIETOOLS 2021b; in the repository since 12/31/2021)
% DJ, 10/29/2026: Update to new pie_struct structure;
% MMP, 11/19/2026: Release comment pass (comment-only).
% DJ, 11/23/2026: Compare the spatial variables in any order;
% MMP, 11/27/2026: copvar, cdopvar renamed opvar, dopvar, with the
%   functions of the rename table; names only.


% List of the operators to compare.
PI_list = {'T', 'Tw', 'Tu',...
           'A', 'Bw', 'Bu',...
           'Cz','Cy', 'Dzw',...
           'Dzu','Dyw', 'Dyu',...
           'K', 'Kw', 'Ku'};                                                % DJ, 10/29/2026

logval = 0;

if isa(objA,'pie_struct')&&isa(objB,'pie_struct')
    % The same variables and domains, in any order
    [is_var,loc] = ismember(objA.space.vars,objB.space.vars);
    if numel(objA.space.vars)~=numel(objB.space.vars) || ~all(is_var) ||...
            ~isequal(objA.dom,objB.dom(loc,:))                              % DJ, 11/23/2026
        logval = 1;
    end
    for i=PI_list
        opA = objA.(i{:});      opB = objB.(i{:});
        if size(opA,1)==0 && size(opB,1)==0
            % Neither of the operators has been specified, or both have no rows.
            continue
        end
        % Operators of different types (e.g. opvar and double), or acting
        % on different spaces, are not equal.
        try is_equal = logical(opA==opB);
        catch, is_equal = false;
        end
        if ~all(is_equal)
            logval = 1;
        end
    end
else
    logval = 1;
end

logval = ~(logval);
end