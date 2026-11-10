function int_order = get_int_order(intmat)
% INT_ORDER = GET_INT_ORDER(INTMAT) returns the order in which the
% variables must be integrated in the term of a functional specified by
% INTMAT.
%
% INPUTS
% - intmat:     1 x 2*d array of integers between 0 and d+1, specifying the
%               lower and upper boundaries of the intervals over which to
%               integrate d variables, as for a row of the field intmat of
%               an 'intop' object;
%
% OUTPUTS
% - int_order:  1 x d array of integers specifying for each of the
%               variables in what order it is integrated. Here,
%               int_order(k) = q implies that variable k must be
%               integrated qth;
%
% EXAMPLE
% Consider d=3, with the two terms specified by the rows of
%   intmat = [0,2,0,4,2,4; 0,4,0,1,0,2]
% specifying the integrals
%   int_{[a,b]^3} chi_{[a,s2]}(s1) chi_{[a,b]}(s2) chi_{[s2,b]} ... ds
%    + int_{[a,b]^3} chi_{[a,b]}(s1) chi_{[a,s1]}(s2) chi_{[a,s2]} ... ds
% representing
%   int_{a}^{b} int_{a}^{s2} int_{s2}^{b} ... ds3 ds1 ds2
%    + int_{a}^{b} int_{a}^{s1} int_{a}^{s2} ... ds3 ds2 ds1
% In the first term, variable s2 must be integrated last, as it is
% integrated over the full domain [a,b]. Variables s1 and s3 can be
% integrated in either order. In the second term, variable s3 must be
% integrated first, then s2, and then s1. Calling this function for the
% first and for the second row, we get
%   int_order = [1,2,1]   and   int_order = [3,2,1],
% respectively.

% DJ, 10/15/2026: Initial coding;

% Split the intervals into lower and upper boundaries
d = numel(intmat)/2;
Aj = intmat(1:2:end);
Bj = intmat(2:2:end);
% Declare the order of the variables in term j
lvls_j = zeros(1,d);
full_int_vars = [0,d+1];   % variables of which a full integral is taken
var_rem = 1:d;             % variables of which interval of integration is still variable
lvl = 1;                    % order assigned to variables over which full integral is taken
while ~isempty(var_rem)
    % Find all variables that are integrated over the full domain
    % at level "lvl"
    is_full_int = ismember(Aj,full_int_vars) & ismember(Bj,full_int_vars);
    new_vars = var_rem(is_full_int);
    lvls_j(new_vars) = lvl;
    % Add the fully integrated variables to the list
    full_int_vars = [full_int_vars,new_vars];
    var_rem = setdiff(var_rem,new_vars);
    Aj = Aj(~is_full_int);      Bj = Bj(~is_full_int);
    % Move to the next level of integrals
    lvl = lvl+1;
    if lvl>d+1
        error("Order of integrals appears to be ambiguous, something is wrong....")
    end
end
% The order of integration should be inner integral first
int_order = max(lvls_j)-lvls_j+1;

end
