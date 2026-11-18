function logval = eq(objA, objB)
% logval = eq(objA,objB) or objA==objB tests whether two "npie_struct"
% objects represent the same PIE.
%
% INPUT
% - objA, objB:    npie_struct class objects defining PIEs.
%
% OUTPUT
% - logval: true if the PIEs involve the same spatial variables and domain,
%           and each of the operators is the same, false otherwise.
%
% NOTES
% Identical operators are equal. Otherwise, operators are compared by
% taking their difference, and checking whether it is zero. Operators of
% different types, or acting on different function spaces, are not equal.
% For polynomials (polyopvar), the difference may not be recognized as zero
% (e.g. for nonlinear terms), in which case the PIEs are not considered
% equal.

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
% DJ, 10/29/2026: Initial coding;

PI_list = {'T','Tw','Tu',...
           'f','Bw','Bu',...
           'gz','Dzw','Dzu',...
           'gy','Dyw','Dyu',...
           'K','Kw','Ku'};

if ~isa(objA,'npie_struct') || ~isa(objB,'npie_struct')
    logval = false;
    return
end
logval = isequal(objA.space.vars,objB.space.vars) && isequal(objA.dom,objB.dom);
for name = PI_list
    opA = objA.(name{1});
    opB = objB.(name{1});
    if isnumeric(opA) && isempty(opA) && isnumeric(opB) && isempty(opB)
        % Neither of the operators has been specified.
        continue
    elseif isequal(opA,opB)
        continue
    elseif ~strcmp(class(opA),class(opB))
        logval = false;
    elseif isa(opA,'polyopvar') || isa(opA,'tensopmat')
        % The difference is only defined for operators on the same spaces.
        try
            if isa(opA,'polyopvar') && size(opA,1)==0 && size(opB,1)==0
                is_zero = true;     % no outputs
            elseif isa(opA,'polyopvar')
                [~,is_zero] = clean_poly(opA-opB);
            else
                [~,is_zero] = clean_op(opA-opB);
            end
        catch
            is_zero = false;
        end
        logval = logval && is_zero;
    else
        logval = logval && isequal(opA,opB);
    end
    if ~logval
        return
    end
end

end
