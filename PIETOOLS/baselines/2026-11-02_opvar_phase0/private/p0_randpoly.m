function x = p0_randpoly(rows, q, vnames)
% X = P0_RANDPOLY(ROWS,Q,VNAMES) returns a ROWS x Q 'polynomial' with
% standard normal coefficients on the monomials of total degree <= 2 in the
% variables named in VNAMES, or a ROWS x Q 'double' when VNAMES is empty.
% The stream is the caller's (P0_APPLY seeds it).
if isempty(vnames)
    x = randn(rows, q);
    return
end
v = p0_vars(vnames);
nv = numel(vnames);
mons = polynomial(1);
for a = 1:nv,   mons = [mons; v(a)];    end                                 %#ok<AGROW>
for a = 1:nv
    for b = a:nv,   mons = [mons; v(a)*v(b)];   end                         %#ok<AGROW>
end
cols = cell(1,q);
for c = 1:q
    cols{c} = randn(rows, numel(mons)) * mons;
end
x = [cols{:}];
end
