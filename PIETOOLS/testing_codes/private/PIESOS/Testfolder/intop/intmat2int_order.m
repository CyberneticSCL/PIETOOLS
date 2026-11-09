function int_order = intmat2int_order(intmat)
% INT_ORDER = INTMAT2INT_ORDER(INTMAT) returns, for each spatial direction
% and each term of an 'intop' object, the order in which the integration
% variables must be integrated, given the field INTMAT of the object.
%
% INPUTS
% - intmat:     1 x N cell, where element i is a qi x 2*di array of
%               integers specifying the intervals of integration along
%               direction i, as for the field intmat of an 'intop' object;
%
% OUTPUTS
% - int_order:  1 x N cell, where element i is a qi x di array of integers.
%               Here, int_order{i}(j,k) = q implies that variable k of
%               direction i in term j must be integrated qth;
%
% EXAMPLE
% Consider N=1 and d=3, with
%   intmat{1} = [0,2,0,4,2,4; 0,4,0,1,0,2]
% specifying the integrals
%   int_{[a,b]^3} chi_{[a,s2]}(s1) chi_{[a,b]}(s2) chi_{[s2,b]} ... ds
%    + int_{[a,b]^3} chi_{[a,b]}(s1) chi_{[a,s1]}(s2) chi_{[a,s2]} ... ds
% representing
%   int_{a}^{b} int_{a}^{s2} int_{s2}^{b} ... ds2 ds1 ds3
%    + int_{a}^{b} int_{a}^{s1} int_{a}^{s2} ... ds3 ds2 ds1
% In the first term, variable s2 must be integrated last, as it is
% integrated over the full domain [a,b]. Variables s1 and s3 can be
% integrated in either order. In the second term, variables s3 must be
% integrated first, then s2, and then s1. We get
%   int_order{1} = [1,2,1; 3,2,1];

% DJ (CC), 10/15/2026: Initial coding;

% Extract the intervals of integration along each spatial direction
N = numel(intmat);
int_order = cell(1,N);
for i=1:N
    % Consider the intervals for the ith variable
    intmat_i = intmat{i};
    [nterms,di] = size(intmat_i);
    di = di/2;
    % 
    var_order_i = zeros(nterms,di);
    for j=1:nterms
        % Split the intervals into lower and upper boundaries
        Aj = intmat_i(j,1:2:end);
        Bj = intmat_i(j,2:2:end);
        % Declare the order of the variables in term j
        lvls_j = zeros(1,di);
        full_int_vars = [0,di+1];   % variables of which a full integral is taken
        var_rem = 1:di;             % varaibles of which interval of integration is still variable
        lvl = 1;                    % order assigned to variables over which full integral is taken
        while ~isempty(var_rem)
            % Find all variables that are integrated over the full domain
            % at level "ord"
            is_full_int = ismember(Aj,full_int_vars) & ismember(Bj,full_int_vars);
            new_vars = var_rem(is_full_int);
            lvls_j(new_vars) = lvl;
            % Add the fully integrated variables to the list
            full_int_vars = [full_int_vars,new_vars];
            var_rem = setdiff(var_rem,new_vars);
            Aj = Aj(~is_full_int);      Bj = Bj(~is_full_int);
            % Move to the next level of integrals
            lvl = lvl+1;
            if lvl>di+1
                error("Order of integrals appears to be ambiguous, something is wrong....")
            end
        end
        % The order of integration should be inner integral first
        var_order_i(j,:) = max(lvls_j)-lvls_j+1;
    end
    int_order{i} = var_order_i;
end

end