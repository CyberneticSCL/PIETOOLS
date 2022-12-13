function logval = isvalidTerm(coeff,termtype,eqntype,stateTerm,...
                stateArg,derOrder,highestDer,intvar,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isvalidTerm.m     PIETOOLS 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function tests whether the term coeff to be added to equation of
% eqntype is a valid operation
% if logval is 0, no issues
% logval is 2, exp has theta when not allowed
% logval is 3, exp has integral with mismatched variable in state argument
% and integration variable
% logval is 4, exp has s or theta when not allowed being added to finite
% dim equations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022  M. Peet, S. Shivakumar, D. Jagt
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
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications



logval=0; % assume its valid, change if the following set of conditions are not satisfied

pvar s theta; 
coeff = convertCharsToStrings(coeff); % convert to string

% if equation is infinite dimensional, then term can be a function of s
if strcmp(eqntype,'infinite')
    if strcmp(termtype,'Multiplier') % if multiplier coeff can only have s, statearg cannot be theta
        coeff = ['(',coeff,')*',stateTerm,'(',stateArg,')'];
        if contains(coeff, 'theta')
            logval = 2;
        end
    elseif strcmp(termtype, 'Integral') %coeff, statearg can have theta if intvar is theta, a and b can be s
        coeff = ['(',coeff,')*',stateTerm,'(',stateArg,')'];
        if contains(coeff, 'theta')&& ~strcmp(intvar, 'd\theta')
            logval = 2;
        end
        if contains(stateArg,'theta')&&~strcmp(intvar, 'd\theta') %x(theta) ds
            logval =3;
        elseif contains(stateArg,'s')&&~strcmp(intvar, 'ds') % x(s) dtheta
            logval =3;
        end
    end
elseif strcmp(eqntype, 'finite')
    if strcmp(termtype,'Multiplier') % if multiplier coeff can only have numbers, statearg cannot be s or theta
        coeff = ['(',coeff,')*',stateTerm,'(',stateArg,')'];
        if contains(coeff, 'theta')|| contains(coeff, 's')
            logval = 4;
        end
    elseif strcmp(termtype,'Integral')
        if (contains(coeff,'s')||contains(stateArg,'s'))&&strcmp(intvar,'d\theta')
            logval=2;
        end
        if (contains(coeff,'theta')||contains(stateArg,'\theta'))&&strcmp(intvar,'ds')
            logval=2;
        end
        if contains(stateArg,'0')||contains(stateArg,'1')
            logval=2;
        end
    else 
        logval=2;
    end
end

end