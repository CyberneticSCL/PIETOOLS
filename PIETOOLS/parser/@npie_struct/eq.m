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
% The spatial variables may be listed in any order. The polynomials f, gz
% and gy are compared after the operators; their difference may not be
% recognized as zero, e.g. for the factors of a nonlinear term in another
% order (x*x_s and x_s*x). The PIEs are then not considered equal, and the
% warning 'npie_struct:eqUndecided' is issued.

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
% DJ, 10/29/2026: Initial coding;

% The operators, and then the polynomials, compared, by field name
PI_list = {'T','Tw','Tu','Bw','Bu','Dzw','Dzu','Dyw','Dyu','K','Kw','Ku',...
           'f','gz','gy'};

% Objects other than 'npie_struct' are not equal; then the variables and
% domain, and each operator in turn
if ~isa(objA,'npie_struct') || ~isa(objB,'npie_struct')
    logval = false;
    return
end
% The same variables and domains, in any order
[is_var,loc] = ismember(objA.space.vars,objB.space.vars);
logval = numel(objA.space.vars)==numel(objB.space.vars) && all(is_var) &&...
            isequal(objA.dom,objB.dom(loc,:));
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
        % A difference of polynomials not shown to be zero does not
        % prove the PIEs differ
        if ~is_zero && isa(opA,'polyopvar')
            warning('npie_struct:eqUndecided',['The polynomials "',name{1},'" of the PIEs could not be shown to be equal;',...
                ' they may still be equal, e.g. with the factors of a nonlinear term in another order (x*x_s and x_s*x).',...
                ' The PIEs are considered not equal.'])
        end
    else
        % Other classes: equal only if identical
        logval = logval && isequal(opA,opB);
    end
    % Stop at the first operator that differs
    if ~logval
        return
    end
end

end
