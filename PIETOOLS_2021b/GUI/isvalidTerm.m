function logval = isvalidTerm(coeff,termtype,eqntype,stateTerm,...
                stateArg,derOrder,highestDer,intvar,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isvalidTerm.m     PIETOOLS 2021b
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