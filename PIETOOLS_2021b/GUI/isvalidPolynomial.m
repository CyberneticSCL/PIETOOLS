function logval = isvalidPolynomial(exp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% isvalidPolynomial.m     PIETOOLS 2021b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function tests whether a mathematical expression exp is a valid
% polynomial in s and theta
% if logval is 0, then exp is a valid polynomial
% logval is 2, exp has unidentified polynomial variable different from s or theta
% logval is 3, exp has polynomial variable with non-integer power
% logval is 4, exp has other error such as mismatched brackets, incorrect
% operands etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEVELOPER LOGS:
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
pvar s theta;

exp = convertCharsToStrings(exp); % convert to string

logval=0;

try
    tmp=eval(exp);
catch err
    if strcmp(err.message,"B must be a constant.")
        logval=3;    
    elseif contains(err.message,"Unrecognized function or variable")
        logval=2;
    elseif strcmp(err.message,"Invalid use of power for polynomials")
        logval=3;
    else
        logval=4;
    end
end
end