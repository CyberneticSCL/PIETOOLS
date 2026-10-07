% Test script for kernel_degree(...)

% How to see monomial degrees of parameters in a dopvar object.
% T = ndopvar2dopvar(Top)
% R1_degmat = full(T.R.R1.degmat)
% R2_degmat = full(T.R.R2.degmat)

% If sos1 is defined as [prog, sos1] = SOS_DP(prog, V_deg, V_mon, x, dom),
% then V_mon = kernel_degree(V_mon) should always be true!

nz = size(V_low.degmat,1); % number of operators.
deg_opts = zeros(1,nz); % max degrees of each operator.

for i=1:nz
    % size = (number of unique terms in polynomial, number of variables in polynomial)
    % row entries = number of times each variable appears in each unique term
    degmat = full(V_low.C.ops{i}.params.degmat);
    sor = sum(degmat,2); % sum over rows = degree of each term.
    deg_opts(i) = max(sor); % max. of sum over rows = degree of polynomial.
end

% deg_opts

% this formula converts max. polynomial degree
[max_deg, arg_deg] = max(deg_opts);
d = V_low.degmat(arg_deg); % degree of FDP in linear form - should always be even for SOS and LFs.
d = ceil(d/2); % degree of FDP in quadratic form - used in formula.

opdeg = (max_deg-1) / (2*d); % this could be negative if deg_opts = 0
opdeg = max(0, opdeg); % maximum polynomial degree.
% could be fractional in case of non-symmetric LF, should be exact for
% SOS_DP.
opdeg = ceil(opdeg)